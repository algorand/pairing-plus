/*!
 This module defines a hash_to_curve trait.
*/

use bls12_381::{ClearH, IsogenyMap, OSSWUMap};
use hash_to_field::{FromRO, HashToField};
use CurveProjective;

type CoordT<PtT> = <PtT as CurveProjective>::Base;

/// Random oracle and injective maps to curve
pub trait HashToCurve {
    /// Random oracle
    fn hash_to_curve<B: AsRef<[u8]>>(msg: B, ciphersuite: u8) -> Self;

    /// Injective encoding
    fn encode_to_curve<B: AsRef<[u8]>>(msg: B, ciphersuite: u8) -> Self;
}

impl<PtT> HashToCurve for PtT
where
    PtT: ClearH + IsogenyMap + OSSWUMap,
    CoordT<PtT>: FromRO,
{
    fn hash_to_curve<B: AsRef<[u8]>>(msg: B, ciphersuite: u8) -> PtT {
        let mut p = {
            let h2f = HashToField::<CoordT<PtT>>::new(msg, Some(&[ciphersuite]));
            let mut tmp = PtT::osswu_map(&h2f.with_ctr(0));
            tmp.add_assign(&PtT::osswu_map(&h2f.with_ctr(1)));
            tmp
        };
        p.isogeny_map();
        p.clear_h();
        p
    }

    fn encode_to_curve<B: AsRef<[u8]>>(msg: B, ciphersuite: u8) -> PtT {
        let mut p = {
            let h2f = HashToField::<CoordT<PtT>>::new(msg, Some(&[ciphersuite]));
            PtT::osswu_map(&h2f.with_ctr(2))
        };
        p.isogeny_map();
        p.clear_h();
        p
    }
}
