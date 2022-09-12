use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Lines};
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
    let mut seqmap: HashMap<String, Vec<usize>> = HashMap::new();
    let filename = args.secondaries.unwrap();

    for l in reader(&filename).lines() {
      let line = l.unwrap();
      let (lineno, qname) = get_qname_and_lineno_secondaries_file(&line);
      match seqmap.get_mut(&qname) {
        Some(vec) => {
          vec.push(lineno);
        }
        None => {
          seqmap.insert(qname.to_string(), vec![lineno]);
        }
      }
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
              match seqmap.get(&qname) {
                Some(list) => {
                  for lineno in list {
                    let (seq, qual) = get_seq_and_qual(&line);
                    println!("{}\t{}\t{}", lineno, seq, qual);
                  }
                }
                None => {}
              };
            }
          }
        }
        Err(e) => eprintln!("{}", e),
      }
    }
  } else if args.pass3 {
    let filename = args.secondaries.unwrap();
    let mut iter = reader(&filename).lines();
    let mut tup = get_next(&mut iter);
    let mut lineno: i64 = 0;
    for l in stdin.lock().lines() {
      match l {
        Ok(line) => {
          if &line[0..1] != "@" {
            if lineno == tup.0 {
              let split = line.split("\t");
              let mut vec: Vec<&str> = split.collect();
              vec[9] = &tup.1;
              vec[10] = &tup.2;
              println!("{}", vec.join("\t"));
              tup = get_next(&mut iter);
            } else {
              println!("{}", line)
            }
            lineno += 1;
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

fn get_next(iter: &mut Lines<Box<dyn BufRead>>) -> (i64, String, String) {
  match iter.next() {
    Some(p) => {
      let str = p.unwrap();
      let mut split = str.split("\t");
      let lineno = split.next().unwrap().parse::<i64>().unwrap();
      let seq = split.next().unwrap();
      let qual = split.next().unwrap();
      (lineno, seq.to_string(), qual.to_string())
    }
    None => (-1, "".to_string(), "".to_string()),
  }
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
