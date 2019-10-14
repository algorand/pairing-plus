/*!
 This module implements hash_to_field and related hashing primitives
 for use with BLS signatures.
*/

use ff::Field;
use hkdf::Hkdf;
use sha2::digest::generic_array::{ArrayLength, GenericArray};
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
        HashToField::from_hash(Sha256::digest(msg.as_ref()), dst)
    }

    /// Create a new struct given a message already hashed with SHA-256
    pub fn from_hash(msg_hash: GenericArray<u8, <Sha256 as Digest>::OutputSize>, dst: Option<&[u8]>) -> HashToField<T> {
        HashToField::<T> {
            msg_hashed: Hkdf::<Sha256>::extract(dst, msg_hash.as_ref()).0,
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
