/*!
 This module defines a hash_to_curve trait.
*/

use bls12_381::{ClearH, IsogenyMap, OSSWUMap};
use hash_to_field::{hash_to_field, ExpandMsg, FromRO};
use map_to_curve::MapToCurve;
use CurveProjective;
use SubgroupCheck;

type CoordT<PtT> = <PtT as CurveProjective>::Base;

/// Random oracle and injective maps to curve
pub trait HashToCurve<X>
where
    X: ExpandMsg,
{
    /// Random oracle
    fn hash_to_curve<Mt: AsRef<[u8]>, Dt: AsRef<[u8]>>(msg: Mt, dst: Dt) -> Self;

    /// Injective encoding
    fn encode_to_curve<Mt: AsRef<[u8]>, Dt: AsRef<[u8]>>(msg: Mt, dst: Dt) -> Self;
}

impl<PtT, X> HashToCurve<X> for PtT
where
    PtT: ClearH + IsogenyMap + OSSWUMap,
    <PtT as CurveProjective>::Affine: SubgroupCheck,
    CoordT<PtT>: FromRO,
    X: ExpandMsg,
{
    fn hash_to_curve<Mt: AsRef<[u8]>, Dt: AsRef<[u8]>>(msg: Mt, dst: Dt) -> PtT {
        let u = hash_to_field::<CoordT<PtT>, X>(msg.as_ref(), dst.as_ref(), 2);
        <PtT as MapToCurve::<PtT>>::map2_to_curve(&u[0], &u[1])
    }

    fn encode_to_curve<Mt: AsRef<[u8]>, Dt: AsRef<[u8]>>(msg: Mt, dst: Dt) -> PtT {
        let u = hash_to_field::<CoordT<PtT>, X>(msg.as_ref(), dst.as_ref(), 1);
        <PtT as MapToCurve::<PtT>>::map_to_curve(&u[0])
    }
}
