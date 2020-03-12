/*!
 This module implements hash_to_field and related hashing primitives
 for use with BLS signatures.
*/

use ff::Field;
use hkdf::Hkdf;
use sha2::digest::generic_array::{typenum::Unsigned, ArrayLength, GenericArray};
use sha2::digest::{BlockInput, ExtendableOutput, Input};
use sha2::{Digest, Sha256};
use std::marker::PhantomData;

/// A struct that handles hashing a message to one or more values of T.
#[derive(Debug)]
pub struct HashToField<T> {
    msg_hashed: GenericArray<u8, <Sha256 as Digest>::OutputSize>,
    ctr: u8,
    phantom: PhantomData<T>,
}

impl<T: FromRO> HashToField<T> {
    /// Create a new struct given a message and ciphersuite.
    pub fn new<B: AsRef<[u8]>>(msg: B, dst: Option<&[u8]>) -> HashToField<T> {
        // Unfortunately, Hkdf requires a &[u8]; not obvious to me how to avoid copying.
        // Would be nicer if Hkdf operated on an Iterator over u8 instead...
        let mut msg_padded = Vec::<u8>::with_capacity(msg.as_ref().len() + 1);
        msg_padded.extend_from_slice(msg.as_ref());
        msg_padded.push(0u8);
        HashToField::from_padded(&msg_padded[..], dst)
    }

    /// Create a new struct given a message already padded with a single 0x00 byte
    pub fn from_padded<B: AsRef<[u8]>>(msg_padded: B, dst: Option<&[u8]>) -> HashToField<T> {
        HashToField::<T> {
            msg_hashed: Hkdf::<Sha256>::extract(dst, msg_padded.as_ref()).0,
            ctr: 0,
            phantom: PhantomData::<T>,
        }
    }

    /// Compute the output of the random oracle specified by `ctr`.
    pub fn with_ctr(&self, ctr: u8) -> T {
        T::from_ro(&self.msg_hashed, ctr)
    }
}

/// Iterator that outputs the sequence of field elements corresponding to increasing `ctr` values.
impl<T: FromRO> Iterator for HashToField<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        if self.ctr == 255 {
            None
        } else {
            self.ctr += 1;
            Some(T::from_ro(self.msg_hashed.as_slice(), self.ctr - 1))
        }
    }
}

/// Trait for types implementing expand_message interface for hash_to_field
pub trait ExpandMsg {
    fn expand_message(msg: &[u8], dst: &[u8], len_in_bytes: usize) -> Vec<u8>;
}

/// Placeholder type for implementing expand_message_xof based on a hash function
#[derive(Debug)]
pub struct ExpandMsgXof<HashT> {
    phantom: PhantomData<HashT>,
}

/// ExpandMsgXof implements expand_message_xof for the ExpandMsg trait
impl<HashT> ExpandMsg for ExpandMsgXof<HashT>
where
    HashT: Default + ExtendableOutput + Input,
{
    fn expand_message(msg: &[u8], dst: &[u8], len_in_bytes: usize) -> Vec<u8> {
        HashT::default()
            .chain(msg)
            .chain([
                (len_in_bytes >> 8) as u8,
                len_in_bytes as u8,
                dst.len() as u8,
            ])
            .chain(dst)
            .vec_result(len_in_bytes)
    }
}

/// Placeholder type for implementing expand_message_xmd based on a hash function
#[derive(Debug)]
pub struct ExpandMsgXmd<HashT> {
    phantom: PhantomData<HashT>,
}

/// ExpandMsgXmd implements expand_message_xmd for the ExpandMsg trait
impl<HashT> ExpandMsg for ExpandMsgXmd<HashT>
where
    HashT: Digest + BlockInput,
{
    fn expand_message(msg: &[u8], dst: &[u8], len_in_bytes: usize) -> Vec<u8> {
        let b_in_bytes = <HashT as Digest>::OutputSize::to_usize();
        let ell = (len_in_bytes + b_in_bytes - 1) / b_in_bytes;
        if ell > 255 {
            panic!("ell was too big in expand_message_xmd");
        }
        let b_0 = HashT::new()
            .chain(GenericArray::<u8, <HashT as BlockInput>::BlockSize>::default())
            .chain(msg)
            .chain([
                (len_in_bytes >> 8) as u8,
                len_in_bytes as u8,
                0u8,
                dst.len() as u8,
            ])
            .chain(dst)
            .result();

        let mut b_vals = Vec::<u8>::with_capacity(ell * b_in_bytes);
        // b_1
        b_vals.extend_from_slice(
            HashT::new()
                .chain(&b_0[..])
                .chain([1u8, dst.len() as u8])
                .chain(dst)
                .result()
                .as_ref(),
        );

        for idx in 1..ell {
            // b_0 XOR b_(idx - 1)
            let mut tmp = GenericArray::<u8, <HashT as Digest>::OutputSize>::default();
            b_0.iter()
                .zip(&b_vals[(idx - 1) * b_in_bytes..idx * b_in_bytes])
                .enumerate()
                .for_each(|(jdx, (b0val, bi1val))| tmp[jdx] = b0val ^ bi1val);
            b_vals.extend_from_slice(
                HashT::new()
                    .chain(tmp)
                    .chain([(idx + 1) as u8, dst.len() as u8])
                    .chain(dst)
                    .result()
                    .as_ref(),
            );
        }

        b_vals.truncate(len_in_bytes);
        b_vals
    }
}

/// Trait implementing hashing to a field or extension.
pub trait FromRO {
    /// from_ro gives the result of hash_to_field(msg, ctr) when input = H(msg).
    fn from_ro<B: AsRef<[u8]>>(input: B, ctr: u8) -> Self;
}

/// Generic implementation for non-extension fields having a BaseFromRO impl.
impl<T: BaseFromRO> FromRO for T {
    fn from_ro<B: AsRef<[u8]>>(input: B, ctr: u8) -> T {
        T::base_from_ro(input.as_ref(), ctr, 1)
    }
}

/// Implements the loop body of hash_to_base from hash-to-curve draft.
pub trait BaseFromRO: Field {
    /// The length of the HKDF output used to hash to a field element.
    type Length: ArrayLength<u8>;

    /// Convert piece of HKDF output to field element
    fn from_okm(okm: &GenericArray<u8, <Self as BaseFromRO>::Length>) -> Self;

    /// Returns the value from the inner loop of hash_to_field by
    /// hashing twice, calling sha_to_base on each, and combining the result.
    fn base_from_ro(msg_hashed: &[u8], ctr: u8, idx: u8) -> Self {
        let mut result = GenericArray::<u8, <Self as BaseFromRO>::Length>::default();
        let h = Hkdf::<Sha256>::from_prk(msg_hashed).unwrap();
        // "H2C" || I2OSP(ctr, 1) || I2OSP(idx, 1)
        let info = [72, 50, 67, ctr, idx];
        h.expand(&info, &mut result).unwrap();
        <Self as BaseFromRO>::from_okm(&result)
    }
}
