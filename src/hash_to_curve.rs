/*!
 This module defines a hash_to_curve trait.
*/

use bls12_381::{ClearH, IsogenyMap, OSSWUMap};
use hash_to_field::{hash_to_field, ExpandMsg, FromRO};
use CurveProjective;

type CoordT<PtT> = <PtT as CurveProjective>::Base;

/// Random oracle and injective maps to curve
pub trait HashToCurve {
    /// Random oracle
    fn hash_to_curve<X: ExpandMsg>(msg: &[u8], dst: &[u8]) -> Self;

    /// Injective encoding
    fn encode_to_curve<X: ExpandMsg>(msg: &[u8], dst: &[u8]) -> Self;
}

impl<PtT> HashToCurve for PtT
where
    PtT: ClearH + IsogenyMap + OSSWUMap,
    CoordT<PtT>: FromRO,
{
    fn hash_to_curve<X: ExpandMsg>(msg: &[u8], dst: &[u8]) -> PtT {
        let mut p = {
            let u = hash_to_field::<CoordT<PtT>, X>(msg, dst, 2);
            let mut tmp = PtT::osswu_map(&u[0]);
            tmp.add_assign(&PtT::osswu_map(&u[1]));
            tmp
        };
        p.isogeny_map();
        p.clear_h();
        p
    }

    fn encode_to_curve<X: ExpandMsg>(msg: &[u8], dst: &[u8]) -> PtT {
        let mut p = {
            let u = hash_to_field::<CoordT<PtT>, X>(msg, dst, 1);
            PtT::osswu_map(&u[0])
        };
        p.isogeny_map();
        p.clear_h();
        p
    }
}
