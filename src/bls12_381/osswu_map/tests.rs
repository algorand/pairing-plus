/*!
Tests for osswu map
*/

use super::OSSWUMap;
use bls12_381::{Fq, Fq2, FqRepr, G1, G2};
use ff::{Field, PrimeField};
use rand_core::SeedableRng;
//use rand::{thread_rng, Rand};

/// check that the point (X : Y : Z)==(X/Z^2, Y/Z^3) is on E: y^2 = x^3 + ELLP_A * x + ELLP_B
fn check_g_prime<F: Field>(x: &F, y: &F, z: &F, a: &F, b: &F) {
    let lhs = {
        // y^2
        let mut tmp = *y;
        tmp.square();
        tmp
    };

    let rhs = {
        // x^3 + A x z^4 + B z^6
        let mut zsq = *z;
        zsq.square();

        let mut z4 = zsq;
        z4.square();

        let mut tmp1 = *x;
        tmp1.square();
        tmp1.mul_assign(x); // x^3

        let mut tmp2 = *x;
        tmp2.mul_assign(&z4);
        tmp2.mul_assign(a);
        tmp1.add_assign(&tmp2); // + A x z^4

        tmp2 = z4;
        tmp2.mul_assign(&zsq);
        tmp2.mul_assign(b);
        tmp1.add_assign(&tmp2); // + B z^6

        tmp1
    };

    assert_eq!(lhs, rhs);
}

fn check_g1_prime(x: &Fq, y: &Fq, z: &Fq) {
    use super::g1::{ELLP_A, ELLP_B};
    check_g_prime(x, y, z, &ELLP_A, &ELLP_B);
}

fn check_g2_prime(x: &Fq2, y: &Fq2, z: &Fq2) {
    use super::g2::{ELLP_A, ELLP_B};
    check_g_prime(x, y, z, &ELLP_A, &ELLP_B);
}

#[test]
fn test_osswu_g1() {
    // exceptional case: zero
    let p = G1::osswu_map(&Fq::zero());
    let G1 { x, y, z } = &p;
    let xo = Fq::from_repr(FqRepr([
        0x6144f0e146df0250u64,
        0x9e9fd4264a7edcbau64,
        0x519289c2e473a9c7u64,
        0xfc9e9c179c1c484fu64,
        0x1bde5cc11dc20ba5u64,
        0x119d96b86f8b3b8bu64,
    ]))
    .unwrap();
    let yo = Fq::from_repr(FqRepr([
        0x2c26d31ff8057aa2u64,
        0x9f824897b954500eu64,
        0xd6b1bcf4165f3575u64,
        0x8d267d9b89fb2b31u64,
        0x905bde90d4b39d8au64,
        0x8327183f6473933u64,
    ]))
    .unwrap();
    let zo = Fq::from_repr(FqRepr([
        0xfe7db859f2cb453fu64,
        0x8e55cb15e9aab878u64,
        0x51fe89284e4d926au64,
        0x9a148b96ab3e6941u64,
        0xa3857e1ea7b2289du64,
        0xdf088f08f205e3u64,
    ]))
    .unwrap();
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g1_prime(x, y, z);

    // exceptional case: sqrt(-1/XI) (positive)
    let excp = Fq::from_repr(FqRepr([
        0x7cc51062bde821b8u64,
        0x88b69520ee5c57fbu64,
        0x46edbdd403fc310u64,
        0x12f01df4948d09ffu64,
        0xdb38f4a9a3d71bdau64,
        0x1f7462c8b6cbf74u64,
    ]))
    .unwrap();
    let p = G1::osswu_map(&excp);
    let G1 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g1_prime(x, y, z);

    // exceptional case: sqrt(-1/XI) (negative)
    let excp = Fq::from_repr(FqRepr([
        0x3d39ef9d421788f3u64,
        0x95f56addc2f7a804u64,
        0x62c1f6c3b6713313u64,
        0x51872d905ef808c0u64,
        0x6fe2b30c9f7490fdu64,
        0x1809cbbdae132725u64,
    ]))
    .unwrap();
    let p = G1::osswu_map(&excp);
    let G1 { x, y, z } = &p;
    let myo = {
        let mut tmp = yo;
        tmp.negate();
        tmp
    };
    assert_eq!(x, &xo);
    assert_eq!(y, &myo);
    assert_eq!(z, &zo);
    check_g1_prime(x, y, z);

    let u = Fq::from_repr(FqRepr([
        0xd4e2aa3bbf9a8255u64,
        0xa79f2ece3390978cu64,
        0x48c1a8fdff541ebau64,
        0x2b17303f8af1ec82u64,
        0x86657cd3fc3d08b5u64,
        0x14f05da1ad4eddc8u64,
    ]))
    .unwrap();
    let xo = Fq::from_repr(FqRepr([
        0xb8e5b32b10dd26f7u64,
        0x8a114aa4ef26ad27u64,
        0xad97709b49ae7c62u64,
        0x9bc765ec50b53945u64,
        0xae99d020a70ca4feu64,
        0x1803cbf9bd2e3815u64,
    ]))
    .unwrap();
    let yo = Fq::from_repr(FqRepr([
        0x498ec4b38b052163u64,
        0xdfb4b3c21c64a917u64,
        0xa6ad223eeba44938u64,
        0xa564373b4a3b1d49u64,
        0x4f3ba7671555ba8eu64,
        0x141f3b7a3a3bc9a1u64,
    ]))
    .unwrap();
    let zo = Fq::from_repr(FqRepr([
        0xc75f9dc8b69d09eeu64,
        0x80824ef4608083ceu64,
        0xfcd339725e80194au64,
        0xda50cf8999450757u64,
        0x35da50fd75b53f96u64,
        0xade87be1822999bu64,
    ]))
    .unwrap();
    let p = G1::osswu_map(&u);
    let G1 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g1_prime(x, y, z);

    let u = Fq::from_repr(FqRepr([
        0xdfad7422a0bab889u64,
        0x4a70b9f85b2c6f5au64,
        0xc042f72ce88d22f5u64,
        0x5be4f1d4b77bef62u64,
        0x99207c0238d7ab04u64,
        0x6135a609e9aad26u64,
    ]))
    .unwrap();
    let xo = Fq::from_repr(FqRepr([
        0xc43f22e4c5179aa6u64,
        0x90750edf071b3149u64,
        0xddd1fb0b077b1269u64,
        0xf5cef22203523563u64,
        0x6c65968a7d59fffcu64,
        0x9ced6809e9858aeu64,
    ]))
    .unwrap();
    let yo = Fq::from_repr(FqRepr([
        0xdc8a4684944f05adu64,
        0x413c96b605fa42c3u64,
        0x9a8be21d1e63b20eu64,
        0xa86a30f86322e838u64,
        0x77ac5472b64d30abu64,
        0x53119aa51ffec68u64,
    ]))
    .unwrap();
    let zo = Fq::from_repr(FqRepr([
        0xa36fa20f6ddcdbfdu64,
        0x517e8ce7336e879au64,
        0xba98cb9cd4519e1eu64,
        0x7537ed7e920203a5u64,
        0xab59f2690f27e4d9u64,
        0x14fac872814de6e3u64,
    ]))
    .unwrap();
    let p = G1::osswu_map(&u);
    let G1 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g1_prime(x, y, z);

    let u = Fq::from_repr(FqRepr([
        0xaf50b546edfc358au64,
        0x3f1897a2f38a122eu64,
        0xdad7bf8fa9eb51beu64,
        0x34c9f03ed6c4ba66u64,
        0x9ee6db517906e388u64,
        0x1097781715e5c672u64,
    ]))
    .unwrap();
    let xo = Fq::from_repr(FqRepr([
        0x8f0c1b27b7d153a1u64,
        0xef591e984e7736c9u64,
        0x7eb7353e36c7a10eu64,
        0xa13c0d70a7f3a5a0u64,
        0x84e37fc496ea7683u64,
        0xfe619171ecfcbd6u64,
    ]))
    .unwrap();
    let yo = Fq::from_repr(FqRepr([
        0xdc73edc70e39fe42u64,
        0x62ecf6761ff9930fu64,
        0x18187783faab9a4cu64,
        0xaf8e80ddfe379c09u64,
        0xfc8ef86e038520aau64,
        0xc5fca1550de691du64,
    ]))
    .unwrap();
    let zo = Fq::from_repr(FqRepr([
        0x5b66b6ee03f15298u64,
        0x89237edcc40aed57u64,
        0x37259c742eca1bb1u64,
        0xe70fee0572e60397u64,
        0x22fce25b7e2597b9u64,
        0x18e223a3b11df7a4u64,
    ]))
    .unwrap();
    let p = G1::osswu_map(&u);
    let G1 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g1_prime(x, y, z);

    let u = Fq::from_repr(FqRepr([
        0xea84b00658419fc4u64,
        0xdc23cabb1c5bedd0u64,
        0x51b2c9560f33a8d5u64,
        0xdce76c736ec4a3d3u64,
        0xaed02316b6641449u64,
        0x17c2c631ba5d8bebu64,
    ]))
    .unwrap();
    let xo = Fq::from_repr(FqRepr([
        0x4387a325ed54b1d1u64,
        0x9e27b0edabd4fe91u64,
        0xca40b0c21fecd54u64,
        0x7fb2ac0251eee168u64,
        0x89a3fb041cc9ad83u64,
        0x163ba2f38efc6de4u64,
    ]))
    .unwrap();
    let yo = Fq::from_repr(FqRepr([
        0x8e4021829edb5acau64,
        0x69f8104daa66ea3eu64,
        0x8e0604fc190b8ad0u64,
        0x661e41dc536ff246u64,
        0x1838aaa49432898fu64,
        0x899f70d9a4252d5u64,
    ]))
    .unwrap();
    let zo = Fq::from_repr(FqRepr([
        0x2e88745aee5b0da3u64,
        0x5ce92018233a731fu64,
        0x2fac5fa03579f6f7u64,
        0x69c2227c1dbcf7b4u64,
        0x65aded420fb38ca4u64,
        0x24327b6cd1e6b84u64,
    ]))
    .unwrap();
    let p = G1::osswu_map(&u);
    let G1 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g1_prime(x, y, z);

    let mut rng = rand_xorshift::XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    for _ in 0..32 {
        let input = Fq::random(&mut rng);
        let p = G1::osswu_map(&input);
        let G1 { x, y, z } = &p;
        check_g1_prime(x, y, z);
    }
}

#[test]
fn test_osswu_g2() {
    let c0 = Fq::from_repr(FqRepr([0xb1e40u64, 0x0u64, 0x0u64, 0x0u64, 0x0u64, 0x0u64])).unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xb9fefffffffbf5ebu64,
        0x1eabfffeb153ffffu64,
        0x6730d2a0f6b0f624u64,
        0x64774b84f38512bfu64,
        0x4b1ba7b6434bacd7u64,
        0x1a0111ea397fe69au64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x28d6c99ea4383807u64,
        0x59cc5836c91ef30fu64,
        0xa87d216900801408u64,
        0x2610ff4c3c3f9eb1u64,
        0x4f4b3ea32be995fcu64,
        0xdc6721ebe6be37u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x60e0254afc5ba93bu64,
        0x407c6124b57df4cu64,
        0xf8f3c1f44b0f8c7au64,
        0xb96a3df0badd28fau64,
        0x3d04c58bb5e6260u64,
        0x12b21ca35569a3eeu64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([0xf0u64, 0x0u64, 0x0u64, 0x0u64, 0x0u64, 0x0u64])).unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xb9feffffffffa8cbu64,
        0x1eabfffeb153ffffu64,
        0x6730d2a0f6b0f624u64,
        0x64774b84f38512bfu64,
        0x4b1ba7b6434bacd7u64,
        0x1a0111ea397fe69au64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&Fq2::zero());
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g2_prime(x, y, z);

    let c0 = Fq::from_repr(FqRepr([0x76980u64, 0x0u64, 0x0u64, 0x0u64, 0x0u64, 0x0u64])).unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x3b4c00u64,
        0x0u64,
        0x0u64,
        0x0u64,
        0x0u64,
        0x0u64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0xe24baa0a898b47e0u64,
        0x92afb1b88e09c84cu64,
        0xf16d677192b7b78au64,
        0xab1dd12189c47c0eu64,
        0xc30f74ce786d38e9u64,
        0xcc49de633f05c98u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x936dda4aedcab1e1u64,
        0x8261a18f1038bdbu64,
        0xc08dea79dde085du64,
        0x9002d76a3ed1ffd2u64,
        0x185ab763985ff885u64,
        0xbab7cc25639665u64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([0x2d0u64, 0x0u64, 0x0u64, 0x0u64, 0x0u64, 0x0u64])).unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xb9feffffffffa9bbu64,
        0x1eabfffeb153ffffu64,
        0x6730d2a0f6b0f624u64,
        0x64774b84f38512bfu64,
        0x4b1ba7b6434bacd7u64,
        0x1a0111ea397fe69au64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&Fq2::one());
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g2_prime(x, y, z);

    let m1 = {
        let mut tmp = Fq2::one();
        tmp.negate();
        tmp
    };
    let p = G2::osswu_map(&m1);
    let myo = {
        let mut tmp = yo;
        tmp.negate();
        tmp
    };
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &myo);
    assert_eq!(z, &zo);
    check_g2_prime(x, y, z);

    let c0 = Fq::from_repr(FqRepr([
        0xd4e2aa3bbf9a8255u64,
        0xa79f2ece3390978cu64,
        0x48c1a8fdff541ebau64,
        0x2b17303f8af1ec82u64,
        0x86657cd3fc3d08b5u64,
        0x14f05da1ad4eddc8u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x472f6df27fe7c94du64,
        0xea72d4e6f4f06693u64,
        0xd1a89c5e84e6d193u64,
        0xab80a6a3842df525u64,
        0x46e112ac0a450ea4u64,
        0x171441a6d04ca8a9u64,
    ]))
    .unwrap();
    let u = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x19399d7e6e728efau64,
        0x9223ea49b3a6685bu64,
        0xb0535eeb3e0be8eeu64,
        0xccdd7c2ed7a70c2du64,
        0x192ab8f31b9bb432u64,
        0xc0b207783a7fe8au64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xc65c4431a6496c30u64,
        0x8542454973283f10u64,
        0xa7808bb40eebf6b9u64,
        0x683e0aad6e74a5a0u64,
        0x2076b05de214ef02u64,
        0xe039ae7c29d2022u64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0xcb28446e62179d9eu64,
        0xa280a992df73998eu64,
        0x2d5291422919d305u64,
        0x418c865e205bc0c6u64,
        0xf8d1e5e8c38550acu64,
        0xee2df0d5e07448fu64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xaa7c2684fe2fcc6eu64,
        0x99a983385cb3106fu64,
        0x37ad3280cb8a1519u64,
        0x5a4308b2de7f901du64,
        0xf2f74d4b44fadc7cu64,
        0x6ac1c85e32f4edcu64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x9e31b14df8456862u64,
        0xb09d54057305d0eau64,
        0x7d4ec28cf63bbd66u64,
        0x1817c2139c736f55u64,
        0x7fd9f027c2ed4347u64,
        0x18d33c46e9efe1f7u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x4da85b1219f0aa69u64,
        0x9eb5f7883c8356b6u64,
        0x9d27373105a8522fu64,
        0x5be18ff40be45f19u64,
        0x9b693bc483f0f59fu64,
        0x922c5bef1fc118cu64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&u);
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g2_prime(x, y, z);

    let c0 = Fq::from_repr(FqRepr([
        0xdfad7422a0bab889u64,
        0x4a70b9f85b2c6f5au64,
        0xc042f72ce88d22f5u64,
        0x5be4f1d4b77bef62u64,
        0x99207c0238d7ab04u64,
        0x6135a609e9aad26u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x34f124763e7deb00u64,
        0xa285e8e52a9cf5f5u64,
        0x3463f5943127700cu64,
        0xeea0ef2a7244c951u64,
        0xeeedf7205412c6a4u64,
        0x3ac7d4da624f424u64,
    ]))
    .unwrap();
    let u = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x470a0eb4c6ea41dau64,
        0x38fc102a7ac96c4bu64,
        0xf12cc75f43f16fau64,
        0x1ae7110401d2bf60u64,
        0xabcdd7ccae9a680au64,
        0x7a6102bf5d97c9cu64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x184b0324bbf4ec25u64,
        0x14e6a614c88543ebu64,
        0x11b6dadcb855c02eu64,
        0x45d1bc1a7b21bf38u64,
        0x6e9811b7292cbe35u64,
        0x20c43c3e504b49du64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x6c3bf81fbed884beu64,
        0xc9913eda07951808u64,
        0x74fe400d891f0d2fu64,
        0x66459536249908u64,
        0xc6dd9a1dd87e2749u64,
        0x15ab62367cef7a16u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xf0f9b256ffd9ffd0u64,
        0xc481edd39780ca8fu64,
        0x5ea12e0601bcb0adu64,
        0x92ffe49990fd9032u64,
        0xacc33c14b83593a7u64,
        0x5048b4e608e7595u64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x85e0d5489736e2b4u64,
        0x6c5118e2091d88f0u64,
        0x8b41f404e6916df1u64,
        0xda99a9546f39acf9u64,
        0x57587e3b4ed7340du64,
        0x170ef6f0827380fcu64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xdd61a360bf21c990u64,
        0xe87c9a8fbef8edfeu64,
        0x674f970b3d82e9b8u64,
        0xb3f831e1eabbf03bu64,
        0xcee9367de3ca318u64,
        0x160a61c5ad6a3ff3u64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&u);
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g2_prime(x, y, z);

    let c0 = Fq::from_repr(FqRepr([
        0xaf50b546edfc358au64,
        0x3f1897a2f38a122eu64,
        0xdad7bf8fa9eb51beu64,
        0x34c9f03ed6c4ba66u64,
        0x9ee6db517906e388u64,
        0x1097781715e5c672u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x9c0ae60f939506a8u64,
        0xa4ef9b76946849beu64,
        0x2d7708869060ff0cu64,
        0xbd6d915e7952a21du64,
        0xbfa926b829513c7eu64,
        0x1732337eace2d016u64,
    ]))
    .unwrap();
    let u = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x80330decf209c0f9u64,
        0x9c3c443d2148943cu64,
        0x7b012833fbb8d302u64,
        0xc46b5c5bdffaf903u64,
        0xdc32da48bd881df2u64,
        0xf7a0d745e96ee8cu64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x1f77ac75a53cb01du64,
        0x331ccd087fe7e20u64,
        0xc798a6624c5c2657u64,
        0x318fdef5c6a03aaeu64,
        0x75d649c08a4329b5u64,
        0xd8461734f2b818du64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x979dcacafe864046u64,
        0x58d5a8feda1b9c68u64,
        0x2dc7fe64d491eb67u64,
        0x1e4eda823cf7dd15u64,
        0x307aa36319482608u64,
        0x2251c003acfa5u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x71eae2c9a83ace48u64,
        0x94e7e90db7f89c6eu64,
        0xeb1f1e2094f14b12u64,
        0x4c44debec0dd26f4u64,
        0x78fd720e9efd3821u64,
        0x145c52e57606ffa5u64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x9b1eaa292b0b8d6cu64,
        0xf3556f782b80156au64,
        0x7232a60dfcf45578u64,
        0xda283bc794f1c552u64,
        0x72e449993919e49au64,
        0xdd03753cbb62029u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x5e1300be109addd0u64,
        0xf9a438110153ac6fu64,
        0x3f16da21234b7dfeu64,
        0x668a29f291c491ccu64,
        0xb007536e7f23b656u64,
        0x1435472d4037af40u64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&u);
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g2_prime(x, y, z);

    let c0 = Fq::from_repr(FqRepr([
        0xea84b00658419fc4u64,
        0xdc23cabb1c5bedd0u64,
        0x51b2c9560f33a8d5u64,
        0xdce76c736ec4a3d3u64,
        0xaed02316b6641449u64,
        0x17c2c631ba5d8bebu64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xb2577499ede5f632u64,
        0xca3d6ab753b878fu64,
        0x1833b9b48c4d08cdu64,
        0x9df66243f1e33375u64,
        0xeecbfb9b9c09d227u64,
        0x7a4a6b660e99b12u64,
    ]))
    .unwrap();
    let u = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x9b719651c4c746e6u64,
        0xbd438453f89d2adcu64,
        0x22116768f501742eu64,
        0x51174b39ab6bc2cu64,
        0xe1c665b1e5c63de6u64,
        0x1842adaf28baae5u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x6b54949d6f96dbcfu64,
        0xa915298df9efc27au64,
        0x3439428ca0b987e5u64,
        0x61ea03ec041d8965u64,
        0x86c6f8125dc0bbc2u64,
        0xddb31de92a06828u64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x4b0f0747efbfaa2bu64,
        0xf5501541912d865bu64,
        0x977c499198af07aau64,
        0xd446b9a7fad8f3b9u64,
        0x4badb10ce5e47e33u64,
        0x907aa22e7410129u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xdd46de8193ef06f6u64,
        0x4cce40b2f67a7123u64,
        0xd843215fe615c9b4u64,
        0x8d839820983c4b41u64,
        0x16042cc81ac4edddu64,
        0x7eff2aeb396734au64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x23f96d2a1601cb7u64,
        0x6a074b0a3175cbfcu64,
        0x28a4ab30815e16a1u64,
        0x1030979d8436dd2eu64,
        0xb43ad04879add9d4u64,
        0x522b59175626baau64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x6992705ff971d0dau64,
        0x295c53f6b1faaa69u64,
        0xe07009934bc1022eu64,
        0x47e2a110d26f261u64,
        0x1721f26639694182u64,
        0x15dba187573a86c3u64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&u);
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);
    check_g2_prime(x, y, z);

    let mut rng = rand_xorshift::XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    for _ in 0..32 {
        let input = Fq2::random(&mut rng);
        let p = G2::osswu_map(&input);
        let G2 { x, y, z } = &p;
        check_g2_prime(x, y, z);
    }
}
