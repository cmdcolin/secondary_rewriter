use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
  #[clap(long, required(false), takes_value(false))]
  pass1: bool,

  #[clap(long, required(false), takes_value(false))]
  pass2: bool,

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

  if args.pass1 {
    for l in stdin.lock().lines() {
      match l {
        Ok(line) => {
          if &line[0..1] != "@" {
            let flag = get_flags(&line);
            if flag & 256 > 0 && flag & 2048 == 0 {
              println!("{}", line)
            }
          }
        }
        Err(e) => eprintln!("{}", e),
      }
    }
  } else if args.pass2 {
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

fn get_seq_and_qual(s: &str) -> (String, String) {
  let mut iter = s.match_indices('\t');
  let i = match iter.nth(8) {
    Some(v) => v.0 + 1,
    None => panic!("Error at line {}", s),
  };
  let l = match iter.nth(0) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  };
  let seq = &s[i..l];

  let i2 = l + 1;
  let l2 = match iter.nth(0) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  };
  let qual = &s[i2..l2];

  (String::from(seq), String::from(qual))
}

fn rewrite_seq_qual(s: &str, seq: &str, qual: &str) -> String {
  let mut iter = s.match_indices('\t');
  let i = match iter.nth(8) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  };

  let j = match iter.nth(1) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  };
  let start = &s[0..i];
  let end = &s[j..];
  format!("{}\t{}\t{}{}", start, seq, qual, end)
}

fn get_flags(s: &str) -> u16 {
  let mut iter = s.match_indices('\t');
  let j = match iter.nth(0) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  };
  let k = match iter.nth(0) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  };
  let flags = &s[j + 1..k];
  flags.parse::<u16>().unwrap()
}

fn get_qname_and_flags(s: &str) -> (String, u16) {
  let mut iter = s.match_indices('\t');
  let j = match iter.nth(0) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  };
  let qname = &s[0..j];
  let k = match iter.nth(0) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  };
  let flags = &s[j + 1..k];
  let f = flags.parse::<u16>().unwrap();
  (String::from(qname), f)
}

fn get_qname(s: &str) -> String {
  let mut iter = s.match_indices('\t');
  let j = match iter.nth(0) {
    Some(v) => v.0,
    None => panic!("Error at line {}", s),
  };
  String::from(&s[0..j])
}
