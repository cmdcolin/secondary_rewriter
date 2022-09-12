use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::error::Error;
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
}

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

fn main() -> Result<(), Box<dyn Error>> {
  let stdin = io::stdin();
  let args = Args::parse();

  // if multiple secondaries with same QNAME, store each in VEC
  let mut seqmap: HashMap<String, Vec<String>> = HashMap::new();
  let filename = args.secondaries.unwrap();

  for l in reader(&filename).lines() {
    let line = l.unwrap();
    let qname = String::from(get_qname(&line));
    match seqmap.get_mut(&qname) {
      Some(vec) => {
        vec.push(line);
      }
      None => {
        seqmap.insert(qname.to_string(), vec![line]);
      }
    }
  }

  for l in stdin.lock().lines() {
    match l {
      Ok(line) => {
        if &line[0..1] != "@" {
          let (qname, flag) = get_qname_and_flags(&line);

          // if it passes this check, it is a primary alignment which we can get SEQ+QUAL from.
          // reason: it's not secondary because input to this program uses -F256 to
          // filter out secondaries and not supplementary, which do have SEQ but can be
          // hard clipped)
          if flag & 2048 == 0 && seqmap.contains_key(&qname) {
            let (seq, qual) = get_seq_and_qual(&line);

            // output secondaries (this implies out-of-order lines now, requires the
            // pipe-to-sort
            match seqmap.get(&qname) {
              Some(list) => {
                for secondary_line in list {
                  let mut vec: Vec<&str> = secondary_line.split("\t").collect();
                  vec[9] = &seq;
                  vec[10] = &qual;
                  println!("{}", vec.join("\t"));
                }
              }
              None => {}
            };
          }
        }
        // output primary and header lines as usual
        println!("{}", line)
      }
      Err(e) => eprintln!("{}", e),
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

fn get_seq_and_qual(s: &str) -> (&str, &str) {
  let mut iter = s.match_indices('\t');
  let i = match_pos(&mut iter, 8, s) + 1;
  let l = match_pos(&mut iter, 0, s);
  let seq = &s[i..l];

  let i2 = l + 1;
  let l2 = match_pos(&mut iter, 0, s);
  let qual = &s[i2..l2];

  (seq, qual)
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

fn get_qname(s: &str) -> &str {
  let mut iter = s.match_indices('\t');
  let j = match_pos(&mut iter, 0, s);
  &s[0..j]
}
