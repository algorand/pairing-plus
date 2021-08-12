/*!
 This module defines the `MapToCurve` trait and implements it
 in terms of `ClearH + IsogenyMap + OSSWUMap`.
*/

use bls12_381::{ClearH, IsogenyMap, OSSWUMap};
use CurveProjective;

pub trait MapToCurve<PtT>
where
	PtT: CurveProjective,
{
	fn map_to_curve(p1: &PtT::Base) -> PtT;
	fn map2_to_curve(p1: &PtT::Base, p2: &PtT::Base) -> PtT;
}

impl<PtT> MapToCurve<PtT> for PtT
where
    PtT: ClearH + IsogenyMap + OSSWUMap,
{
	fn map_to_curve(p1: &PtT::Base) -> PtT {
        let mut p = PtT::osswu_map(p1);
        p.isogeny_map();
        p.clear_h();
        p
	}

	fn map2_to_curve(p1: &PtT::Base, p2: &PtT::Base) -> PtT {
        let mut p = {
            let mut tmp = PtT::osswu_map(p1);
            tmp.add_assign(&PtT::osswu_map(p2));
            tmp
        };
        p.isogeny_map();
        p.clear_h();
        p
	}
}
