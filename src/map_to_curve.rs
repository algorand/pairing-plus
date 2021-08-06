use bls12_381::{ClearH, IsogenyMap, OSSWUMap};
use CurveProjective;

pub trait MapToCurve<PtT>
where
	PtT: CurveProjective,
{
	fn map_one_field_element_to_curve(p1: &PtT::Base) -> PtT;
	fn map_two_field_elements_to_curve(p1: &PtT::Base, p2: &PtT::Base) -> PtT;
}

impl<PtT> MapToCurve<PtT> for PtT
where
    PtT: ClearH + IsogenyMap + OSSWUMap,
{
	fn map_one_field_element_to_curve(p1: &PtT::Base) -> PtT {
        println!("map_one_field_element_to_curve called");
        let mut p = PtT::osswu_map(p1);
        p.isogeny_map();
        p.clear_h();
        p
	}

	fn map_two_field_elements_to_curve(p1: &PtT::Base, p2: &PtT::Base) -> PtT {
        println!("map_two_field_elements_to_curve called");
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