mod g1 {
    use pairing::bls12_381::*;
    use pairing::CurveAffine;
    use pairing::CurveProjective;
    use rand::{Rand, Rng, SeedableRng, XorShiftRng};
    use ff::PrimeField;
    #[bench]
    fn bench_g1_mul_shamir(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G1, G1, Fr, Fr)> = (0..SAMPLES)
            .map(|_| {
                (
                    (G1::rand(&mut rng)),
                    (G1::rand(&mut rng)),
                    Fr::rand(&mut rng),
                    Fr::rand(&mut rng),
                )
            })
            .collect();
        let mut count = 0;
        b.iter(|| {
            let tmp = CurveProjective::mul_shamir(v[count].0, v[count].1, v[count].2, v[count].3);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g1_mul_assign(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G1, Fr)> = (0..SAMPLES)
            .map(|_| (G1::rand(&mut rng), Fr::rand(&mut rng)))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = v[count].0;
            tmp.mul_assign(v[count].1);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g1_membership(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<G1> = (0..SAMPLES).map(|_| G1::rand(&mut rng)).collect();

        let mut count = 0;
        b.iter(|| {
            let tmp = v[count].into_affine().membership_check();
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g1_mul_assign_sec(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        // mul_assign_sec ensures constant time for a same base point
        // and various scalars
        let p = G1::rand(&mut rng);
        let v: Vec<Fr> = (0..SAMPLES).map(|_| Fr::rand(&mut rng)).collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = p.clone();
            tmp.mul_assign_sec(v[count]);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g1_mul_cofactor(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<G1Affine> = (0..SAMPLES)
            .map(|_| (G1::rand(&mut rng)).into_affine())
            .collect();

        let mut count = 0;
        b.iter(|| {
            let tmp = v[count];
            tmp.scale_by_cofactor();
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g1_add_assign(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G1, G1)> = (0..SAMPLES)
            .map(|_| (G1::rand(&mut rng), G1::rand(&mut rng)))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = v[count].0;
            tmp.add_assign(&v[count].1);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g1_add_assign_mixed(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G1, G1Affine)> = (0..SAMPLES)
            .map(|_| (G1::rand(&mut rng), G1::rand(&mut rng).into()))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = v[count].0;
            tmp.add_assign_mixed(&v[count].1);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g1_doubling(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G1, G1Affine)> = (0..SAMPLES)
            .map(|_| (G1::rand(&mut rng), G1::rand(&mut rng).into()))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = v[count].0;
            tmp.double();
            count = (count + 1) % SAMPLES;
            tmp
        });
    }
    #[bench]
    fn bench_cast_string_to_e1(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = rand::thread_rng();
        let mut inputstr: [u8; 48] = [0; 48];
        let mut count = 0;
        b.iter(|| {
            for x in inputstr.iter_mut() {
                *x = rng.gen();
            }

            let tmp = G1Affine::cast_string_to_e1(inputstr);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }
    #[bench]
    fn hash_to_e1(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut v: Vec<String> = Vec::new();
        for _i in 0..SAMPLES {
            let s = rand::thread_rng()
                .gen_ascii_chars()
                .take(10)
                .collect::<String>();
            v.push(s);
        }

        let mut count = 0;
        b.iter(|| {
            let tmp = G1Affine::hash_to_e1(v[count].as_bytes());
            count = (count + 1) % SAMPLES;
            tmp
        });
    }
    #[bench]
    fn hash_to_g1(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut v: Vec<String> = Vec::new();
        for _i in 0..SAMPLES {
            let s = rand::thread_rng()
                .gen_ascii_chars()
                .take(10)
                .collect::<String>();
            v.push(s);
        }

        let mut count = 0;
        b.iter(|| {
            let tmp = G1Affine::hash_to_g1(v[count].as_bytes());
            count = (count + 1) % SAMPLES;
            tmp
        });
    }
    #[bench]
    fn hash_to_g1_const(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut v: Vec<String> = Vec::new();
        for _i in 0..SAMPLES {
            let s = rand::thread_rng()
                .gen_ascii_chars()
                .take(10)
                .collect::<String>();
            v.push(s);
        }

        let mut count = 0;
        b.iter(|| {
            let tmp = G1Affine::hash_to_g1_const(v[count].as_bytes());
            count = (count + 1) % SAMPLES;
            tmp
        });
    }
    #[bench]
    fn solve_for_y(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;
            let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
        let mut v: Vec<Fq> = Vec::new();
        for _i in 0..SAMPLES {
            let s = Fq::rand(&mut rng);
            v.push(s);
        }

        let mut count = 0;
        b.iter(|| {
            let tmp = G1Affine::rhs_g1(&v[count]);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g1_sw_encode(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;
            let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
        let mut v: Vec<Fq> = Vec::new();
        for _i in 0..SAMPLES {
            let s = Fq::rand(&mut rng);
            v.push(s);
        }

        let mut count = 0;
        b.iter(|| {
            let tmp = G1Affine::g1_sw_encode(v[count]);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }
    #[bench]
    fn bench_g1_get_point(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;
            let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
        let mut v: Vec<Fq> = Vec::new();
        for _i in 0..SAMPLES {
            let s = Fq::rand(&mut rng);
            v.push(s);
        }

        let mut count = 0;
        b.iter(|| {
            let tmp = G1Affine::get_point_from_x(v[count], true);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }
    #[bench]
    fn bench_conversion(b: &mut ::test::Bencher) {


        pub const SQRT_NEG_THREE: &str = "1586958781458431025242759403266842894121773480562120986020912974854563298150952611241517463240701";
        pub const SQRT_NEG_THREE_MIN_ONE_DIV_TWO: &str="793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350";

        let mut count = 0;
        b.iter(|| {

                let sqrt_neg_three = Fq::from_str(&SQRT_NEG_THREE).unwrap();
            //count = (count + 1) % SAMPLES;
            //tmp
        });
    }


}

mod g2 {
    use pairing::bls12_381::*;
    use pairing::CurveAffine;
    use pairing::CurveProjective;
    use rand::{Rand, Rng, SeedableRng, XorShiftRng};

    #[bench]
    fn bench_g2_mul_shamir(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G2, G2, Fr, Fr)> = (0..SAMPLES)
            .map(|_| {
                (
                    (G2::rand(&mut rng)),
                    (G2::rand(&mut rng)),
                    Fr::rand(&mut rng),
                    Fr::rand(&mut rng),
                )
            })
            .collect();
        let mut count = 0;
        b.iter(|| {
            let tmp = CurveProjective::mul_shamir(v[count].0, v[count].1, v[count].2, v[count].3);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g2_mul_assign(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G2, Fr)> = (0..SAMPLES)
            .map(|_| (G2::rand(&mut rng), Fr::rand(&mut rng)))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = v[count].0;
            tmp.mul_assign(v[count].1);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g2_membership(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<G2> = (0..SAMPLES).map(|_| G2::rand(&mut rng)).collect();

        let mut count = 0;
        b.iter(|| {
            let tmp = v[count].into_affine().membership_check();
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g2_mul_assign_sec(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
        let p: G2 = G2::rand(&mut rng);
        let v: Vec<Fr> = (0..SAMPLES).map(|_| Fr::rand(&mut rng)).collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = p.clone();
            tmp.mul_assign_sec(v[count]);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g2_mul_cofactor(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<G2Affine> = (0..SAMPLES)
            .map(|_| (G2::rand(&mut rng)).into_affine())
            .collect();

        let mut count = 0;
        b.iter(|| {
            let tmp = v[count];
            tmp.scale_by_cofactor();
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g2_add_assign(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G2, G2)> = (0..SAMPLES)
            .map(|_| (G2::rand(&mut rng), G2::rand(&mut rng)))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = v[count].0;
            tmp.add_assign(&v[count].1);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g2_add_assign_mixed(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G2, G2Affine)> = (0..SAMPLES)
            .map(|_| (G2::rand(&mut rng), G2::rand(&mut rng).into()))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = v[count].0;
            tmp.add_assign_mixed(&v[count].1);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_g2_doubling(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        let v: Vec<(G2, G2Affine)> = (0..SAMPLES)
            .map(|_| (G2::rand(&mut rng), G2::rand(&mut rng).into()))
            .collect();

        let mut count = 0;
        b.iter(|| {
            let mut tmp = v[count].0;
            tmp.double();
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn bench_cast_string_to_e2(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut rng = rand::thread_rng();
        let mut inputstr: [u8; 96] = [0; 96];
        let mut count = 0;
        b.iter(|| {
            for x in inputstr.iter_mut() {
                *x = rng.gen();
            }

            let tmp = G1Affine::cast_string_to_e2(inputstr);
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn hash_to_e2(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut v: Vec<String> = Vec::new();
        for _i in 0..SAMPLES {
            let s = rand::thread_rng()
                .gen_ascii_chars()
                .take(10)
                .collect::<String>();
            v.push(s);
        }

        let mut count = 0;
        b.iter(|| {
            //let ref input_str = &v[count];
            let tmp = G2Affine::hash_to_e2(v[count].as_bytes());
            count = (count + 1) % SAMPLES;
            tmp
        });
    }

    #[bench]
    fn hash_to_g2(b: &mut ::test::Bencher) {
        const SAMPLES: usize = 1000;

        let mut v: Vec<String> = Vec::new();
        for _i in 0..SAMPLES {
            let s = rand::thread_rng()
                .gen_ascii_chars()
                .take(10)
                .collect::<String>();
            v.push(s);
        }

        let mut count = 0;
        b.iter(|| {
            //let ref input_str = &v[count];
            let tmp = G2Affine::hash_to_g2(v[count].as_bytes());
            count = (count + 1) % SAMPLES;
            tmp
        });
    }
}
