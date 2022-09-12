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

  #[clap(long, required(false), takes_value(false))]
  pass1: bool,

  #[clap(long, required(false), takes_value(false))]
  pass2: bool,

  #[clap(long, required(false), takes_value(false))]
  pass3: bool,
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

fn main() -> Result<(), Box<dyn Error>> {
  let stdin = io::stdin();
  let args = Args::parse();

  if args.pass1 {
    let mut lineno = 0;
    for l in stdin.lock().lines() {
      match l {
        Ok(line) => {
          if &line[0..1] != "@" {
            let flag = get_flags(&line);
            if flag & 256 > 0 {
              println!("{}\t{}", lineno, line)
            }
            lineno += 1;
          }
        }
        Err(e) => eprintln!("{}", e),
      }
    }
  } else if args.pass2 {
    let mut seqmap: HashMap<String, Option<(String, String)>> = HashMap::new();
    let filename = args.secondaries.unwrap();
    for l in reader(&filename).lines() {
      let line = l.unwrap();
      let (_lineno, qname) = get_qname_and_lineno_secondaries_file(&line);
      seqmap.insert(qname.to_string(), None);
    }

    for l in stdin.lock().lines() {
      match l {
        Ok(line) => {
          if &line[0..1] != "@" {
            let (qname, flag) = get_qname_and_flags(&line);

            // don't print secondaries again till the end, which are obtained from
            // secondaries.txt. filter out supplementary alignments also: look for
            // actual primary alignment
            if flag & 256 == 0 && flag & 2048 == 0 && seqmap.contains_key(&qname) {
              seqmap.insert(qname, Some(get_seq_and_qual(&line)));
            }
          }
        }
        Err(e) => eprintln!("{}", e),
      }
    }

    for l in reader(&filename).lines() {
      let line = l.unwrap();
      let (_lineno, qname) = get_qname_and_lineno_secondaries_file(&line);
      match seqmap.get(&qname).unwrap() {
        Some(e) => {
          let (seq, qual) = e;
          println!("{}", rewrite_seq_qual(&line, &seq, &qual));
        }
        None => eprintln!("Primary record for QNAME not found: {}", qname),
      }
    }
  } else if args.pass3 {
    let filename = args.secondaries.unwrap();
    let mut seqmap = HashMap::new();
    for l in reader(&filename).lines() {
      let line = l.unwrap();
      let (lineno, rest) = line.split_once("\t").unwrap();
      seqmap.insert(lineno.parse::<usize>().unwrap(), String::from(rest));
    }

    let mut i: usize = 0;
    for l in stdin.lock().lines() {
      match l {
        Ok(line) => {
          if &line[0..1] != "@" {
            if seqmap.contains_key(&i) {
              println!("{}", seqmap.get(&i).unwrap());
            } else {
              println!("{}", line)
            }

            i += 1;
          } else {
            println!("{}", line);
          }
        }
        Err(e) => eprintln!("{}", e),
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
  let i = match_pos(&mut iter, 9, s);
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

// this is from a SAM-like format that has line number of it's position in the original file
// prepended
fn get_qname_and_lineno_secondaries_file(s: &str) -> (usize, String) {
  let mut iter = s.match_indices('\t');
  let j = match_pos(&mut iter, 0, s);
  let k = match_pos(&mut iter, 0, s);
  let lineno = &s[0..j];

  (lineno.parse::<usize>().unwrap(), String::from(&s[j + 1..k]))
}

fn get_flags(s: &str) -> u16 {
  let mut iter = s.match_indices('\t');
  let j = match_pos(&mut iter, 0, s) + 1;
  let k = match_pos(&mut iter, 0, s);
  let flags = &s[j..k];
  flags.parse::<u16>().unwrap()
}
