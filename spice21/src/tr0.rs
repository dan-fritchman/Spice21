use crate::assert::assert;
use crate::spresult::SpResult;
use byteorder::{NativeEndian, ReadBytesExt};
use std::fs::File;
use std::io::prelude::*;
use std::io::SeekFrom;
use std::str;

const VERSION: &[u8] = "9601".as_bytes();
const MAGIC_KEY: &[u8] = "$&%#".as_bytes();
const ENDMAGIC: f32 = 0.99e30; // Generally 1e30, but parsed equality often fails
                               // Sizes of all sorts of tr0 magic-values, all in bytes.
const SIZEOF_INT: usize = 4;
const SIZEOF_FLOAT: usize = 4;
const MAGIC_LEN: usize = 24;
const SIZEOF_BLOCK_SEPARATOR: usize = 20;

/// Independent Variable Type
#[derive(Debug)]
enum IndepType {
    TIME = 1,
    FREQ = 2,
    PARAM = 3,
}
/// Variable datatype for real-values, from DC, tran etc.
#[derive(Debug)]
enum RealVar {
    V = 1,
    I = 8,
}
/// AC Variable Datatypes
#[derive(Debug)]
enum AcVar {
    V_COMPLEX = 1,
    V_DB = 6,
    I_COMPLEX = 8,
    I_DB = 13,
}
/// Datatype Enum
#[derive(Debug)]
enum ValType {
    RE(RealVar),
    AC(AcVar),
}
#[derive(Debug)]
struct Tr0 {
    filename: String,
    header: String,
    indep: i32,           //FIXME: IndepType,
    data_types: Vec<i32>, //FIXME: Vec<ValType>,
    block_sizes: Vec<usize>,
    signal_names: Vec<String>,
    nsig: usize,
    len: usize,
    data_start: usize,
    end: usize,
}
impl Tr0 {
    /// Retrieve values of signal `name`
    fn get(&self, name: &str) -> Option<Vec<f64>> {
        let pos = self.signal_names.iter().position(|x| x == name)?;
        let step = (SIZEOF_FLOAT * (self.nsig - 1)) as i64; // Note `-1` is because read itself moves us 4B
        let end = self.end as u64;
        let mut cursor = (self.data_start + pos * 4) as u64;
        let mut vals: Vec<f64> = vec![];
        let mut block_num = 0;
        let mut next_block_start = (self.data_start + self.block_sizes[0]) as u64;

        // Re-open our file
        let mut file = File::open(&self.filename).unwrap();
        file.seek(SeekFrom::Start(cursor)).unwrap(); 
        // And collect data values 
        loop {
            if cursor >= end {
                break;
            }
            if cursor >= next_block_start {
                cursor = file.seek(SeekFrom::Current(SIZEOF_BLOCK_SEPARATOR as i64)).unwrap();
                next_block_start += (SIZEOF_BLOCK_SEPARATOR + self.block_sizes[block_num]) as u64;
                block_num += 1;
            }
            let v = file.read_f32::<NativeEndian>().unwrap() as f64;
            vals.push(v);
            cursor = file.seek(SeekFrom::Current(step)).unwrap();
        }
        Some(vals)
    }
    fn open(file_name: &str) -> SpResult<Tr0> {
        // Open our file, read its header
        let mut file = File::open(&file_name).unwrap();
        const BUFFER_SIZE: usize = 1024;
        let mut buffer = [0; BUFFER_SIZE];
        // Find the version string
        let _ = file.take(BUFFER_SIZE as u64).read(&mut buffer).unwrap();
        let mut version_start = 0;
        for k in 0..BUFFER_SIZE - VERSION.len() {
            if buffer[k..k + VERSION.len()] == *VERSION {
                version_start = k;
                break;
            }
        }
        assert(version_start).ne(0)?;
        let mut magic_start = 0;
        for k in version_start..BUFFER_SIZE - MAGIC_KEY.len() {
            if buffer[k..k + MAGIC_KEY.len()] == *MAGIC_KEY {
                magic_start = k;
                break;
            }
        }
        assert(magic_start).ne(0)?;
        let block_size_start = magic_start + MAGIC_LEN;
        let data_start = block_size_start + SIZEOF_INT;

        // Convert the header to string, and chop it into signals and datatypes
        let header = str::from_utf8(&buffer[version_start..magic_start]).unwrap();
        let vec: Vec<&str> = header.split("  0  ").collect();
        assert(vec.len()).gt(1)?;
        let sigs = vec[vec.len() - 1];
        let sigs_vec: Vec<&str> = sigs.split_whitespace().collect();
        assert(sigs_vec.len()).gt(1)?;
        assert(sigs_vec.len() % 2).eq(0)?;
        let sig_dtypes: Vec<i32> = sigs_vec[0..sigs_vec.len() / 2].iter().map(|x| x.parse::<i32>().unwrap()).collect();
        let sig_names: Vec<String> = sigs_vec[sigs_vec.len() / 2..sigs_vec.len()].iter().map(|x| x.to_string()).collect();
        assert(sig_dtypes.len()).eq(sigs_vec.len() / 2)?;
        assert(sig_names.len()).eq(sigs_vec.len() / 2)?;

        // Decode the first block-size
        let mut block_size_bytes = &buffer[block_size_start..block_size_start + SIZEOF_INT];
        let mut block_size = block_size_bytes.read_u32::<NativeEndian>().unwrap() as usize;
        let mut block_sizes: Vec<usize> = vec![block_size];

        // Reopen our file (sadly `take` consumed it) to parse the block delineations
        let mut file = File::open(&file_name).unwrap();
        file.seek(SeekFrom::Start((data_start + block_size - SIZEOF_FLOAT) as u64)).unwrap(); 
        loop {
            let last_value = file.read_f32::<NativeEndian>().unwrap();
            // Blocks end by repeating their sizes. (?). Read this in, and check that it matches.
            let past_block_size = file.read_u32::<NativeEndian>().unwrap() as usize;
            assert(past_block_size).eq(block_size)?;
            // Check for the end-of-data magic-value
            if last_value > ENDMAGIC {
                break;
            }
            // Another block present, decode its size. Discard 12B of not-sure-what in between.
            file.seek(SeekFrom::Current((SIZEOF_BLOCK_SEPARATOR - SIZEOF_FLOAT - SIZEOF_INT) as i64))
                .unwrap();
            block_size = file.read_u32::<NativeEndian>().unwrap() as usize;
            block_sizes.push(block_size);
            file.seek(SeekFrom::Current((block_size - SIZEOF_FLOAT) as i64)).unwrap();
        }
        let end = file.seek(SeekFrom::Current(0)).unwrap() as usize - SIZEOF_FLOAT - SIZEOF_INT;
        let nsig = sig_names.len();
        let num_separator_bytes = (block_sizes.len() - 1) * SIZEOF_BLOCK_SEPARATOR;
        // Check we have a valid data-size, modulo the number of signals
        let num_data_vals = end - data_start - num_separator_bytes;
        let mod_ = num_data_vals / SIZEOF_FLOAT % nsig;
        assert(mod_).eq(0)?;
        let len = num_data_vals / SIZEOF_FLOAT / nsig;
        // Loading succeeded. Create and return our struct.
        let tr0 = Tr0 {
            filename: file_name.to_string(),
            header: header.to_string(),
            indep: sig_dtypes[0],
            data_types: sig_dtypes,
            signal_names: sig_names,
            block_sizes,
            data_start,
            nsig,
            end,
            len,
        };
        Ok(tr0)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::spresult::{sperror, TestResult};
    #[test]
    fn test_tr0_read() -> TestResult {
        let file_name = "tbd.tr0";

        let tr0 = Tr0::open(file_name)?;
        // Check all the signals retrieve same-length vectors
        let names = tr0.signal_names.clone();
        for name in names.iter() {
            let vals = tr0.get(name).ok_or(sperror("Signal Not Found"))?;
            assert(vals.len()).eq(tr0.len)?;
        }
        Ok(())
    }
}
