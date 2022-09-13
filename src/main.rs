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
          if flag & 2048 == 0 {
            // output secondaries (this implies out-of-order lines now, requires the
            // pipe-to-sort
            match seqmap.get(qname) {
              Some(list) => {
                let (seq, qual) = get_seq_and_qual(&line);
                let primary_rev = get_flags(&line) & 16;
                for secondary_line in list {
                  let secondary_rev = get_flags(&secondary_line) & 16;
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

fn get_flags(s: &str) -> u16 {
  let mut iter = s.match_indices('\t');
  let j = match_pos(&mut iter, 0, s) + 1;
  let k = match_pos(&mut iter, 0, s);
  let flags = &s[j..k];
  flags.parse::<u16>().unwrap()
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

fn get_qname_and_flags(s: &str) -> (&str, u16) {
  let mut iter = s.match_indices('\t');
  let j = match_pos(&mut iter, 0, s);

  let qname = &s[0..j];
  let k = match_pos(&mut iter, 0, s);
  let flags = &s[j + 1..k];
  let f = flags.parse::<u16>().unwrap();
  (qname, f)
}

fn get_qname(s: &str) -> &str {
  let mut iter = s.match_indices('\t');
  let j = match_pos(&mut iter, 0, s);
  &s[0..j]
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn qname() {
    let q=get_qname("ctgA_3_555_0:0:0_2:0:0_102d\t0\tctgA\t3\t37\t100M\t*\t0\t0\tTTGTTGCGGAGTTGAACAACGGCATTAGGAACACTTCCGTCTCTCACTTTTATACGATTATGATTGGTTCTTTAGCCTTGGTTTAGATTGGTAGTAGTAG\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222\tXT:A:U\tNM:i:0\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:100");

    assert_eq!(q, "ctgA_3_555_0:0:0_2:0:0_102d");
  }

  #[test]
  fn flag() {
    let q=get_flags("ctgA_3_555_0:0:0_2:0:0_102d\t0\tctgA\t3\t37\t100M\t*\t0\t0\tTTGTTGCGGAGTTGAACAACGGCATTAGGAACACTTCCGTCTCTCACTTTTATACGATTATGATTGGTTCTTTAGCCTTGGTTTAGATTGGTAGTAGTAG\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222\tXT:A:U\tNM:i:0\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:100");

    assert_eq!(q, 0);
  }

  #[test]
  fn qname_and_flags() {
    let q=get_qname_and_flags("ctgA_3_555_0:0:0_2:0:0_102d\t0\tctgA\t3\t37\t100M\t*\t0\t0\tTTGTTGCGGAGTTGAACAACGGCATTAGGAACACTTCCGTCTCTCACTTTTATACGATTATGATTGGTTCTTTAGCCTTGGTTTAGATTGGTAGTAGTAG\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222\tXT:A:U\tNM:i:0\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:100");

    assert_eq!(q, ("ctgA_3_555_0:0:0_2:0:0_102d", 0));
  }

  #[test]
  fn seq_and_quals() {
    let q=get_seq_and_qual("ctgA_3_555_0:0:0_2:0:0_102d\t0\tctgA\t3\t37\t100M\t*\t0\t0\tTTGTTGCGGAGTTGAACAACGGCATTAGGAACACTTCCGTCTCTCACTTTTATACGATTATGATTGGTTCTTTAGCCTTGGTTTAGATTGGTAGTAGTAG\t2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222\tXT:A:U\tNM:i:0\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:100");

    assert_eq!(q, ("TTGTTGCGGAGTTGAACAACGGCATTAGGAACACTTCCGTCTCTCACTTTTATACGATTATGATTGGTTCTTTAGCCTTGGTTTAGATTGGTAGTAGTAG", "2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"));
  }
}
