use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
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
    // let mut set = HashSet::new();
    let mut seqmap: HashMap<String, Option<(String, String)>> = HashMap::new();
    // File hosts must exist in current path before this produces output
    //
    //
    let filename = args.secondaries.unwrap();
    if let Ok(lines) = read_lines(&filename) {
      for line in lines {
        let l = line.unwrap();
        let qname = get_qname(&l);
        seqmap.insert(qname.to_string(), None);
      }
    }
    for l in stdin.lock().lines() {
      match l {
        Ok(line) => {
          if &line[0..1] != "@" {
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
    if let Ok(lines) = read_lines(&filename) {
      for line in lines {
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
  }

  Ok(())
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
  P: AsRef<Path>,
{
  let file = File::open(filename)?;
  Ok(io::BufReader::new(file).lines())
}

fn get_seq_and_qual(s: &str) -> (String, String) {
  let mut iter = s.match_indices('\t');
  let i = iter.nth(8).unwrap().0 + 1;
  let l = iter.nth(0).unwrap().0;
  let seq = &s[i..l];

  let i2 = l + 1;
  let l2 = iter.nth(0).unwrap().0;
  let qual = &s[i2..l2];

  (String::from(seq), String::from(qual))
}

fn rewrite_seq_qual(s: &str, seq: &str, qual: &str) -> String {
  let mut iter = s.match_indices('\t');
  let i = iter.nth(8).unwrap().0;
  let j = iter.nth(1).unwrap().0;
  let start = &s[0..i];
  let end = &s[j..];
  format!("{}\t{}\t{}{}", start, seq, qual, end)
}

fn get_flags(s: &str) -> u16 {
  let mut iter = s.match_indices('\t');
  let j = iter.nth(0).unwrap().0;
  let k = iter.nth(0).unwrap().0;
  let flags = &s[j + 1..k];
  flags.parse::<u16>().unwrap()
}

fn get_qname_and_flags(s: &str) -> (String, u16) {
  let mut iter = s.match_indices('\t');
  let j = iter.nth(0).unwrap().0;
  let qname = &s[0..j];
  let k = iter.nth(0).unwrap().0;
  let flags = &s[j + 1..k];
  let f = flags.parse::<u16>().unwrap();
  (String::from(qname), f)
}

fn get_qname(s: &str) -> String {
  let mut iter = s.match_indices('\t');
  let j = iter.nth(0).unwrap().0;
  String::from(&s[0..j])
}
