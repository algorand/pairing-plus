use bls12_381::{self, *};
use std::io::{Error, ErrorKind, Read, Result, Write};
use CurveAffine;
use CurveProjective;
use EncodedPoint;
type Compressed = bool;

/// Serialization support for group elements.
pub trait SerDes: Sized {
    /// Serialize a struct to a writer with a flag of compressness.
    fn serialize<W: Write>(&self, writer: &mut W, compressed: Compressed) -> Result<()>;

    /// Deserialize a struct; also returns a flag
    /// if the element was compressed or not.
    fn deserialize<R: Read>(reader: &mut R) -> Result<(Self, Compressed)>;
}

impl SerDes for G1 {
    /// Convert a G1 point to a blob.
    fn serialize<W: Write>(&self, writer: &mut W, compressed: Compressed) -> Result<()> {
        let t = self.into_affine();

        // convert element into an (un)compressed byte string
        let buf = {
            if compressed {
                let tmp = bls12_381::G1Compressed::from_affine(t);
                tmp.as_ref().to_vec()
            } else {
                let tmp = bls12_381::G1Uncompressed::from_affine(t);
                tmp.as_ref().to_vec()
            }
        };

        // format the output
        writer.write_all(&buf)?;
        Ok(())
    }

    /// Deserialize a G1 element from a blob.
    /// Returns an error if deserialization fails.
    fn deserialize<R: Read>(reader: &mut R) -> Result<(Self, Compressed)> {
        // read into buf of compressed size
        let mut buf = vec![0u8; G1Compressed::size()];
        reader.read_exact(&mut buf)?;

        // check the first bit of buf[0] to decide if the point is compressed
        // or not
        if (buf[0] & 0x80) == 0x80 {
            // first bit is 1 => compressed mode
            // convert the blob into a group element
            let mut g_buf = G1Compressed::empty();
            g_buf.as_mut().copy_from_slice(&buf);
            let g = match g_buf.into_affine() {
                Ok(p) => p.into_projective(),
                Err(e) => return Err(Error::new(ErrorKind::InvalidData, e)),
            };
            Ok((g, true))
        } else if (buf[0] & 0x80) == 0x00 {
            // first bit is 0 => uncompressed mode
            // read the next uncompressed - compressed size
            let mut buf2 = vec![0u8; G1Uncompressed::size() - G1Compressed::size()];
            reader.read_exact(&mut buf2)?;
            // now buf holds the whole uncompressed bytes
            buf.append(&mut buf2);
            // convert the buf into a group element
            let mut g_buf = G1Uncompressed::empty();
            g_buf.as_mut().copy_from_slice(&buf);
            let g = match g_buf.into_affine() {
                Ok(p) => p.into_projective(),
                Err(e) => return Err(Error::new(ErrorKind::InvalidData, e)),
            };
            Ok((g, false))
        } else {
            Err(Error::new(
                ErrorKind::InvalidData,
                "Should never reach here. Something is wrong",
            ))
        }
    }
}

impl SerDes for G2 {
    /// Convert a G2 point to a blob.
    fn serialize<W: Write>(&self, writer: &mut W, compressed: Compressed) -> Result<()> {
        let t = self.into_affine();
        // convert element into an (un)compressed byte string
        let buf = {
            if compressed {
                let tmp = bls12_381::G2Compressed::from_affine(t);
                tmp.as_ref().to_vec()
            } else {
                let tmp = bls12_381::G2Uncompressed::from_affine(t);
                tmp.as_ref().to_vec()
            }
        };

        // format the output
        writer.write_all(&buf)?;
        Ok(())
    }

    /// Deserialize a G2 element from a blob.
    /// Returns an error if deserialization fails.
    fn deserialize<R: Read>(reader: &mut R) -> Result<(Self, Compressed)> {
        // read into buf of compressed size
        let mut buf = vec![0u8; G2Compressed::size()];
        reader.read_exact(&mut buf)?;

        // check the first bit of buf[0] to decide if the point is compressed
        // or not
        if (buf[0] & 0x80) == 0x80 {
            // first bit is 1 => compressed mode
            // convert the buf into a group element
            let mut g_buf = G2Compressed::empty();
            g_buf.as_mut().copy_from_slice(&buf);
            let g = match g_buf.into_affine() {
                Ok(p) => p.into_projective(),
                Err(e) => return Err(Error::new(ErrorKind::InvalidData, e)),
            };
            Ok((g, true))
        } else if (buf[0] & 0x80) == 0x00 {
            // first bit is 0 => uncompressed mode
            // read the next uncompressed - compressed size
            let mut buf2 = vec![0u8; G2Uncompressed::size() - G2Compressed::size()];
            reader.read_exact(&mut buf2)?;
            // now buf holds the whole uncompressed bytes
            buf.append(&mut buf2);
            // convert the buf into a group element
            let mut g_buf = G2Uncompressed::empty();
            g_buf.as_mut().copy_from_slice(&buf);
            let g = match g_buf.into_affine() {
                Ok(p) => p.into_projective(),
                Err(e) => return Err(Error::new(ErrorKind::InvalidData, e)),
            };
            Ok((g, false))
        } else {
            Err(Error::new(
                ErrorKind::InvalidData,
                "Should never reach here. Something is wrong",
            ))
        }
    }
}

#[cfg(test)]
mod serdes_test {
    use super::*;
    #[test]
    fn test_g1_serialization_rand() {
        // G1::zero, compressed
        let g1_zero = G1::zero();
        let mut buf: Vec<u8> = vec![];
        // serialize a G1 element into buffer
        assert!(g1_zero.serialize(&mut buf, true).is_ok());
        assert_eq!(buf.len(), 48, "length of blob is incorrect");
        let (g1_zero_recover, compressed) = G1::deserialize(&mut buf[..].as_ref()).unwrap();
        assert_eq!(compressed, true);
        assert_eq!(g1_zero, g1_zero_recover);

        // G1::one, compressed
        let g1_one = G1::one();
        let mut buf: Vec<u8> = vec![];
        // serialize a G1 element into buffer
        assert!(g1_one.serialize(&mut buf, true).is_ok());
        assert_eq!(buf.len(), 48, "length of blob is incorrect");
        let (g1_one_recover, compressed) = G1::deserialize(&mut buf[..].as_ref()).unwrap();
        assert_eq!(compressed, true);
        assert_eq!(g1_one, g1_one_recover);

        // G1::zero, uncompressed
        let mut buf: Vec<u8> = vec![];
        // serialize a G1 element into buffer
        assert!(g1_zero.serialize(&mut buf, false).is_ok());
        assert_eq!(buf.len(), 96, "length of blob is incorrect");
        let (g1_zero_recover, compressed) = G1::deserialize(&mut buf[..].as_ref()).unwrap();
        assert_eq!(compressed, false);
        assert_eq!(g1_zero, g1_zero_recover);

        // G1::one, uncompressed
        let mut buf: Vec<u8> = vec![];
        // serialize a G1 element into buffer
        assert!(g1_one.serialize(&mut buf, false).is_ok());
        assert_eq!(buf.len(), 96, "length of blob is incorrect");
        let (g1_one_recover, compressed) = G1::deserialize(&mut buf[..].as_ref()).unwrap();
        assert_eq!(compressed, false);
        assert_eq!(g1_one, g1_one_recover);
    }

    #[test]
    fn test_g2_serialization_rand() {
        // G2::zero, compressed
        let g2_zero = G2::zero();
        let mut buf: Vec<u8> = vec![];
        // serialize a G2 element into buffer
        assert!(g2_zero.serialize(&mut buf, true).is_ok());
        assert_eq!(buf.len(), 96, "length of blob is incorrect");
        let (g2_zero_recover, compressed) = G2::deserialize(&mut buf[..].as_ref()).unwrap();
        assert_eq!(compressed, true);
        assert_eq!(g2_zero, g2_zero_recover);

        // G2::one, compressed
        let g2_one = G2::one();
        let mut buf: Vec<u8> = vec![];
        // serialize a G2 element into buffer
        assert!(g2_one.serialize(&mut buf, true).is_ok());
        assert_eq!(buf.len(), 96, "length of blob is incorrect");
        let (g2_one_recover, compressed) = G2::deserialize(&mut buf[..].as_ref()).unwrap();
        assert_eq!(compressed, true);
        assert_eq!(g2_one, g2_one_recover);

        // G2::zero, uncompressed
        let mut buf: Vec<u8> = vec![];
        // serialize a G2 element into buffer
        assert!(g2_zero.serialize(&mut buf, false).is_ok());
        assert_eq!(buf.len(), 192, "length of blob is incorrect");
        let (g2_zero_recover, compressed) = G2::deserialize(&mut buf[..].as_ref()).unwrap();
        assert_eq!(compressed, false);
        assert_eq!(g2_zero, g2_zero_recover);

        // G2::one, uncompressed
        let mut buf: Vec<u8> = vec![];
        // serialize a G2 element into buffer
        assert!(g2_one.serialize(&mut buf, false).is_ok());
        assert_eq!(buf.len(), 192, "length of blob is incorrect");
        let (g2_one_recover, compressed) = G2::deserialize(&mut buf[..].as_ref()).unwrap();

        assert_eq!(compressed, false);
        assert_eq!(g2_one, g2_one_recover);
    }
}
