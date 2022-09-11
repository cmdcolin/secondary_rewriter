use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::str::MatchIndices;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
  #[clap(short, long, required(false), value_parser)]
  secondaries: Option<String>,

  #[clap(long, required(false), takes_value(false))]
  output_only_new_data: bool,
}

// Read normal or compressed files seamlessly
// Uses the presence of a `gz` extension to choose between the two
pub fn reader(filename: &str) -> Box<dyn BufRead> {
  let path = Path::new(filename);
  let file = match File::open(&path) {
    Err(e) => panic!("couldn't open {} {}", path.display(), e),
    Ok(file) => file,
  };

  if path.extension() == Some(OsStr::new("gz")) {
    Box::new(BufReader::with_capacity(128 * 1024, GzDecoder::new(file)))
  } else {
    Box::new(BufReader::with_capacity(128 * 1024, file))
  }
}

fn main() -> Result<(), &'static str> {
  let stdin = io::stdin();
  let args = Args::parse();

  let mut seqmap: HashMap<String, Option<(String, String)>> = HashMap::new();
  let filename = args.secondaries.unwrap();
  for line in reader(&filename).lines() {
    let l = line.unwrap();
    let qname = get_qname(&l);
    seqmap.insert(qname.to_string(), None);
  }
  for l in stdin.lock().lines() {
    match l {
      Ok(line) => {
        if line.len() > 0 && &line[0..1] != "@" {
          let (qname, flag) = get_qname_and_flags(&line);

          // don't print secondaries again till the end, which are obtained from
          // secondaries.txt
          if flag & 256 == 0 {
            // avoid supplementary reads with flag 2048, which may have partial (e.g. hard
            // clipped) sequence. the true primary alignment is !secondary and
            // !supplementary
            if seqmap.contains_key(&qname) && flag & 2048 == 0 {
              seqmap.insert(qname, Some(get_seq_and_qual(&line)));
            }

            // pass thru any other data to stdout as long as the "output_only_new_data" is
            // false
            if !args.output_only_new_data {
              println!("{}", line);
            }
          }
        } else {
          println!("{}", line);
        }
      }
      Err(e) => eprintln!("{}", e),
    }

    for line in reader(&filename).lines() {
      let l = line.unwrap();
      let qname = get_qname(&l);
      match seqmap.get(&qname).unwrap() {
        Some(e) => {
          let (seq, qual) = e;
          println!("{}", rewrite_seq_qual(&l, &seq, &qual));
        }
        None => eprintln!("Primary record for QNAME not found: {}", qname),
      }
    }
  }

  Ok(())
}
#[inline(always)]
fn match_pos(iter: &mut MatchIndices<char>, n: usize, s: &str) -> usize {
  match iter.nth(n) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  }
}

fn get_seq_and_qual(s: &str) -> (String, String) {
  let mut iter = s.match_indices('\t');
  let i = match_pos(&mut iter, 8, s) + 1;
  let l = match_pos(&mut iter, 0, s);
  let seq = &s[i..l];

  let i2 = l + 1;
  let l2 = match_pos(&mut iter, 0, s);
  let qual = &s[i2..l2];

  (String::from(seq), String::from(qual))
}

fn rewrite_seq_qual(s: &str, seq: &str, qual: &str) -> String {
  let mut iter = s.match_indices('\t');
  let i = match_pos(&mut iter, 8, s);
  let j = match_pos(&mut iter, 1, s);

  let start = &s[0..i];
  let end = &s[j..];
  format!("{}\t{}\t{}{}", start, seq, qual, end)
}

fn get_qname_and_flags(s: &str) -> (String, u16) {
  let mut iter = s.match_indices('\t');
  let j = match_pos(&mut iter, 0, s);

  let qname = &s[0..j];
  let k = match_pos(&mut iter, 0, s);
  let flags = &s[j + 1..k];
  let f = flags.parse::<u16>().unwrap();
  (String::from(qname), f)
}

fn get_qname(s: &str) -> String {
  let mut iter = s.match_indices('\t');
  let j = match_pos(&mut iter, 0, s);
  String::from(&s[0..j])
}
