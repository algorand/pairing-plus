// `clippy` is a code linting tool for improving code quality by catching
// common mistakes or strange code patterns. If the `cargo-clippy` feature
// is provided, all compiler warnings are prohibited.
#![cfg_attr(feature = "cargo-clippy", deny(warnings))]
#![cfg_attr(feature = "cargo-clippy", allow(inline_always))]
#![cfg_attr(feature = "cargo-clippy", allow(too_many_arguments))]
#![cfg_attr(feature = "cargo-clippy", allow(unreadable_literal))]
#![cfg_attr(feature = "cargo-clippy", allow(many_single_char_names))]
#![cfg_attr(feature = "cargo-clippy", allow(new_without_default_derive))]
#![cfg_attr(feature = "cargo-clippy", allow(write_literal))]
// Force public structures to implement Debug
#![deny(missing_debug_implementations)]
extern crate byteorder;
//#[macro_use]
extern crate bigint;
extern crate ff;
extern crate rand;
extern crate sha2;

#[cfg(test)]
pub mod tests;

pub mod bls12_381;
mod wnaf;
pub use self::wnaf::Wnaf;

use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr, ScalarEngine, SqrtField};
use std::error::Error;
use std::fmt;

/// An "engine" is a collection of types (fields, elliptic curve groups, etc.)
/// with well-defined relationships. In particular, the G1/G2 curve groups are
/// of prime order `r`, and are equipped with a bilinear pairing function.
pub trait Engine: ScalarEngine {
    /// The projective representation of an element in G1.
    type G1: CurveProjective<
            Engine = Self,
            Base = Self::Fq,
            Scalar = Self::Fr,
            Affine = Self::G1Affine,
        > + From<Self::G1Affine>;

    /// The affine representation of an element in G1.
    type G1Affine: CurveAffine<
            Engine = Self,
            Base = Self::Fq,
            Scalar = Self::Fr,
            Projective = Self::G1,
            Pair = Self::G2Affine,
            PairingResult = Self::Fqk,
        > + From<Self::G1>;

    /// The projective representation of an element in G2.
    type G2: CurveProjective<
            Engine = Self,
            Base = Self::Fqe,
            Scalar = Self::Fr,
            Affine = Self::G2Affine,
        > + From<Self::G2Affine>;

    /// The affine representation of an element in G2.
    type G2Affine: CurveAffine<
            Engine = Self,
            Base = Self::Fqe,
            Scalar = Self::Fr,
            Projective = Self::G2,
            Pair = Self::G1Affine,
            PairingResult = Self::Fqk,
        > + From<Self::G2>;

    /// The base field that hosts G1.
    type Fq: PrimeField + SqrtField;

    /// The extension field that hosts G2.
    type Fqe: SqrtField;

    /// The extension field that hosts the target group of the pairing.
    type Fqk: Field;

    /// Perform a miller loop with some number of (G1, G2) pairs.
    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<
            Item = &'a (
                &'a <Self::G1Affine as CurveAffine>::Prepared,
                &'a <Self::G2Affine as CurveAffine>::Prepared,
            ),
        >;

    /// Perform final exponentiation of the result of a miller loop.
    fn final_exponentiation(&Self::Fqk) -> Option<Self::Fqk>;

    /// Performs a complete pairing operation `(p, q)`.
    fn pairing<G1, G2>(p: G1, q: G2) -> Self::Fqk
    where
        G1: Into<Self::G1Affine>,
        G2: Into<Self::G2Affine>,
    {
        Self::final_exponentiation(&Self::miller_loop(
            [(&(p.into().prepare()), &(q.into().prepare()))].into_iter(),
        ))
        .unwrap()
    }

    fn pairing_product<G1, G2>(p1: G1, q1: G2, p2: G1, q2: G2) -> Self::Fqk
    where
        G1: Into<Self::G1Affine>,
        G2: Into<Self::G2Affine>,
    {
        Self::final_exponentiation(&Self::miller_loop(
            [
                (&(p1.into().prepare()), &(q1.into().prepare())),
                (&(p2.into().prepare()), &(q2.into().prepare())),
            ]
            .into_iter(),
        ))
        .unwrap()
    }

    fn pairing_multi_product(p: &[Self::G1Affine], q: &[Self::G2Affine]) -> Self::Fqk {
        let prep_p: Vec<<Self::G1Affine as CurveAffine>::Prepared> =
            p.iter().map(|v| v.prepare()).collect();
        let prep_q: Vec<<Self::G2Affine as CurveAffine>::Prepared> =
            q.iter().map(|v| v.prepare()).collect();
        let mut pairs = Vec::with_capacity(p.len());
        for i in 0..p.len() {
            pairs.push((&prep_p[i], &prep_q[i]));
        }
        let t = Self::miller_loop(&pairs);
        Self::final_exponentiation(&t).unwrap()
    }
}

/// Projective representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait CurveProjective:
    PartialEq
    + Eq
    + Sized
    + Copy
    + Clone
    + Send
    + Sync
    + fmt::Debug
    + fmt::Display
    + rand::Rand
    + 'static
{
    type Engine: Engine<Fr = Self::Scalar>;
    type Scalar: PrimeField + SqrtField;
    type Base: SqrtField;
    type Affine: CurveAffine<Projective = Self, Scalar = Self::Scalar>;

    /// Returns the additive identity.
    fn zero() -> Self;

    /// Returns a fixed generator of unknown exponent.
    fn one() -> Self;

    /// Determines if this point is the point at infinity.
    fn is_zero(&self) -> bool;

    /// Normalizes a slice of projective elements so that
    /// conversion to affine is cheap.
    fn batch_normalization(v: &mut [Self]);

    /// Checks if the point is already "normalized" so that
    /// cheap affine conversion is possible.
    fn is_normalized(&self) -> bool;

    /// Doubles this element.
    fn double(&mut self);

    /// Adds another element to this element.
    fn add_assign(&mut self, other: &Self);
    /// Adds another element to this element.
    fn add_assign_fake(&mut self, other: &Self);

    /// Subtracts another element from this element.
    fn sub_assign(&mut self, other: &Self) {
        let mut tmp = *other;
        tmp.negate();
        self.add_assign(&tmp);
    }

    /// Adds an affine element to this element.
    fn add_assign_mixed(&mut self, other: &Self::Affine);
    /// Adds an affine element to this element.
    fn add_assign_mixed_fake(&mut self, other: &Self::Affine);

    /// Negates this element.
    fn negate(&mut self);

    /// Performs scalar multiplication of this element.
    fn mul_assign<S: Into<<Self::Scalar as PrimeField>::Repr>>(&mut self, other: S);

    /// Performs scalar multiplication of this element in constant time.
    fn mul_assign_sec<S: Into<<Self::Scalar as PrimeField>::Repr>>(&mut self, other: S);

    /// Converts this element into its affine representation.
    fn into_affine(&self) -> Self::Affine;

    /// Recommends a wNAF window table size given a scalar. Always returns a number
    /// between 2 and 22, inclusive.
    fn recommended_wnaf_for_scalar(scalar: <Self::Scalar as PrimeField>::Repr) -> usize;

    /// Recommends a wNAF window size given the number of scalars you intend to multiply
    /// a base by. Always returns a number between 2 and 22, inclusive.
    fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize;

    /// multiplication with shamir's Trick
    /// compute s1 * p1 + s2 * p2 simultaneously
    fn mul_shamir<S: Into<<Self::Scalar as PrimeField>::Repr>>(
        p1: Self,
        p2: Self,
        s1: S,
        s2: S,
    ) -> Self;
}

/// Affine representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait CurveAffine:
    Copy + Clone + Sized + Send + Sync + fmt::Debug + fmt::Display + PartialEq + Eq + 'static
{
    type Engine: Engine<Fr = Self::Scalar>;
    type Scalar: PrimeField + SqrtField;
    type Base: SqrtField;
    type Projective: CurveProjective<Affine = Self, Scalar = Self::Scalar>;
    type Prepared: Clone + Send + Sync + 'static;
    type Uncompressed: EncodedPoint<Affine = Self>;
    type Compressed: EncodedPoint<Affine = Self>;
    type Pair: CurveAffine<Pair = Self>;
    type PairingResult: Field;

    /// Returns the additive identity.
    fn zero() -> Self;

    /// Returns a fixed generator of unknown exponent.
    fn one() -> Self;

    /// Determines if this point represents the point at infinity; the
    /// additive identity.
    fn is_zero(&self) -> bool;

    /// Negates this element.
    fn negate(&mut self);

    /// Performs scalar multiplication of this element with sliding window.
    fn mul<S: Into<<Self::Scalar as PrimeField>::Repr>>(&self, other: S) -> Self::Projective;

    /// Performs scalar multiplication of this element with mixed addition.
    /// use Montgomery method to achieve constant time
    fn mul_sec<S: Into<<Self::Scalar as PrimeField>::Repr>>(&self, other: S) -> Self::Projective;

    /// Prepares this element for pairing purposes.
    fn prepare(&self) -> Self::Prepared;

    /// Perform a pairing
    fn pairing_with(&self, other: &Self::Pair) -> Self::PairingResult;

    /// Converts this element into its affine representation.
    fn into_projective(&self) -> Self::Projective;

    /// Converts this element into its compressed encoding, so long as it's not
    /// the point at infinity.
    fn into_compressed(&self) -> Self::Compressed {
        <Self::Compressed as EncodedPoint>::from_affine(*self)
    }

    /// Converts this element into its uncompressed encoding, so long as it's not
    /// the point at infinity.
    fn into_uncompressed(&self) -> Self::Uncompressed {
        <Self::Uncompressed as EncodedPoint>::from_affine(*self)
    }

    // this function casts a compressed string into an element on E1
    // it may fail if the compressed string does not correspond to
    // an E1 element
    fn cast_string_to_e1(s: [u8; 48]) -> Result<bls12_381::G1Affine, GroupDecodingError>;

    // this function hashes a string into an E1 element (not gauranteed on G1)
    // the current implementation uses try_and_increment method with
    // cast_string_to_e1 subroutine
    // not a constant time implementation
    fn hash_to_e1(input: &[u8]) -> bls12_381::G1Affine;

    // this function hashes a string into an G1 element
    // ensures constant number of calls
    // not a constant time implementation due to underly math
    fn hash_to_g1_const(input: &[u8]) -> bls12_381::G1;

    // subrouting of hash_to_g1_const
    // uses Shallue-van de Woestijne encoding
    // as described in Fouque-Tibouchi paper
    fn g1_sw_encode(input: bls12_381::Fq) -> bls12_381::G1Affine;

    // given x, compute x^3+b
    fn rhs_g1(x: &bls12_381::Fq) -> bls12_381::Fq;

    // this function hashes a string into a G1 element (gauranteed)
    // uses hash_to_e1 subroutine
    // not a constant time implementation
    fn hash_to_g1(input: &[u8]) -> bls12_381::G1;

    // this function casts a compressed string into an element on E2
    // it may fail if the compressed string does not correspond to
    // an E2 element
    fn cast_string_to_e2(s: [u8; 96]) -> Result<bls12_381::G2Affine, GroupDecodingError>;

    // this function hashes a string into an E2 element (not gauranteed on G2)
    // the current implementation uses try_and_increment method with
    // cast_string_to_e2 subroutine
    // not a constant time implementation
    fn hash_to_e2(input: &[u8]) -> bls12_381::G2Affine;

    // this function hashes a string into a G2 element (gauranteed)
    // uses hash_to_e2 subroutine
    // not a constant time implementation
    fn hash_to_g2(input: &[u8]) -> bls12_381::G2;

    /// multiplication of many points
    /// compute s1 * p1 + ... + sn * pn simultaneously
    fn sum_of_products(bases: &[Self], scalars: &[&[u64; 4]]) -> Self::Projective;

    /// Find the optimal window for running Pippinger's algorithm; preprogrammed values
    fn find_pippinger_window(num_components: usize) -> usize;
    /// Find the optimal window for running Pippinger's algorithm; computed values via an estimate of running time
    fn find_pippinger_window_via_estimate(num_components: usize) -> usize;
    /// multiplication of many points with Pippinger's algorithm of window size w
    /// compute s1 * p1 + ... + sn * pn simultaneously
    fn sum_of_products_pippinger(
        bases: &[Self],
        scalars: &[&[u64; 4]],
        window: usize,
    ) -> Self::Projective;
    /// multiplication of many points with Pippinger's algorithm of window size w
    /// compute s1 * p1 + ... + sn * pn simultaneously
    fn sum_of_products_pippinger_projective_only(
        bases: &[Self],
        scalars: &[&[u64; 4]],
        window: usize,
    ) -> Self::Projective;
    /// multiplication of many points with Pippinger's algorithm of window size w
    /// compute s1 * p1 + ... + sn * pn simultaneously
    fn sum_of_products_pippinger_affine_only(
        bases: &[Self],
        scalars: &[&[u64; 4]],
        window: usize,
    ) -> Self::Projective;
    /// multiplication of many points with Pippinger's algorithm of window size w
    /// compute s1 * p1 + ... + sn * pn simultaneously
    fn sum_of_products_pippinger_fake_math(
        bases: &[Self],
        scalars: &[&[u64; 4]],
        window: usize,
    ) -> Self::Projective;

    /// multiplication of many points with precompuation
    /// compute s1 * p1 + ... + sn * pn simultaneously
    /// assuming  pre[j*256+i] = (\sum_{b such that bth bit of i is 1} 2^{32i}) * bases[j] for each j and i in 0..256
    fn sum_of_products_precomp_256(
        bases: &[Self],
        scalars: &[&[u64; 4]],
        pre: &[Self],
    ) -> Self::Projective;

    /// pre[0] becomes (2^64) * self, pre[1]  becomes (2^128) * self, and pre[2] (becomes 2^196) * self
    fn precomp_3(&self, pre: &mut [Self]);

    /// Performs scalar multiplication of this element,
    /// assuming pre = [(2^64)*self, (2^128)*self, (2^192)*self]
    fn mul_precomp_3<S: Into<<Self::Scalar as PrimeField>::Repr>>(
        &self,
        other: S,
        pre: &[Self],
    ) -> Self::Projective;

    /// pre[i] becomes (\sum_{b such that bth bit of i is 1} 2^{32i}) * self for i in 0..25
    fn precomp_256(&self, pre: &mut [Self]);

    /// Performs scalar multiplication of this element,
    /// assuming  pre[i] = (\sum_{b such that bth bit of i is 1} 2^{32i}) * self for i in 0..256
    fn mul_precomp_256<S: Into<<Self::Scalar as PrimeField>::Repr>>(
        &self,
        other: S,
        pre: &[Self],
    ) -> Self::Projective;
}

/// An encoded elliptic curve point, which should essentially wrap a `[u8; N]`.
pub trait EncodedPoint:
    Sized + Send + Sync + AsRef<[u8]> + AsMut<[u8]> + Clone + Copy + 'static
{
    type Affine: CurveAffine;

    /// Creates an empty representation.
    fn empty() -> Self;

    /// Returns the number of bytes consumed by this representation.
    fn size() -> usize;

    /// Converts an `EncodedPoint` into a `CurveAffine` element,
    /// if the encoding represents a valid element.
    fn into_affine(&self) -> Result<Self::Affine, GroupDecodingError>;

    /// Converts an `EncodedPoint` into a `CurveAffine` element,
    /// without guaranteeing that the encoding represents a valid
    /// element. This is useful when the caller knows the encoding is
    /// valid already.
    ///
    /// If the encoding is invalid, this can break API invariants,
    /// so caution is strongly encouraged.
    fn into_affine_unchecked(&self) -> Result<Self::Affine, GroupDecodingError>;

    /// Creates an `EncodedPoint` from an affine point, as long as the
    /// point is not the point at infinity.
    fn from_affine(affine: Self::Affine) -> Self;
}

/// An error that may occur when trying to decode an `EncodedPoint`.
#[derive(Debug)]
pub enum GroupDecodingError {
    /// The coordinate(s) do not lie on the curve.
    NotOnCurve,
    /// The element is not part of the r-order subgroup.
    NotInSubgroup,
    /// One of the coordinates could not be decoded
    CoordinateDecodingError(&'static str, PrimeFieldDecodingError),
    /// The compression mode of the encoded element was not as expected
    UnexpectedCompressionMode,
    /// The encoding contained bits that should not have been set
    UnexpectedInformation,
}

impl Error for GroupDecodingError {
    fn description(&self) -> &str {
        match *self {
            GroupDecodingError::NotOnCurve => "coordinate(s) do not lie on the curve",
            GroupDecodingError::NotInSubgroup => "the element is not part of an r-order subgroup",
            GroupDecodingError::CoordinateDecodingError(..) => "coordinate(s) could not be decoded",
            GroupDecodingError::UnexpectedCompressionMode => {
                "encoding has unexpected compression mode"
            }
            GroupDecodingError::UnexpectedInformation => "encoding has unexpected information",
        }
    }
}

impl fmt::Display for GroupDecodingError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match *self {
            GroupDecodingError::CoordinateDecodingError(description, ref err) => {
                write!(f, "{} decoding error: {}", description, err)
            }
            _ => write!(f, "{}", self.description()),
        }
    }
}
