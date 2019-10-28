/*!
Tests for osswu map
*/

use super::OSSWUMap;
use bls12_381::{Fq, Fq2, FqRepr, G1, G2};
use ff::{Field, PrimeField};
use rand::{thread_rng, Rand};

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

    let mut rng = thread_rng();
    for _ in 0..32 {
        let input = Fq::rand(&mut rng);
        let p = G1::osswu_map(&input);
        let G1 { x, y, z } = &p;
        check_g1_prime(x, y, z);
    }
}

#[test]
fn test_osswu_g2() {
    let c0 = Fq::from_repr(FqRepr([
        0xb9fefffffff8412bu64,
        0x1eabfffeb153ffffu64,
        0x6730d2a0f6b0f624u64,
        0x64774b84f38512bfu64,
        0x4b1ba7b6434bacd7u64,
        0x1a0111ea397fe69au64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x1bbe3a0d7831b22du64,
        0x57d06f44dc428a47u64,
        0x5f2c8926e75ab5b9u64,
        0x950d9410caa33cf3u64,
        0x76f86a9f629c9333u64,
        0x1109ddf9a05fe2c7u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x42344d91d4b86b64u64,
        0xb6664db1979e52e6u64,
        0x1db9054108eda6f2u64,
        0xf5714ec6406123e6u64,
        0x8adfb4a7ca0f35e2u64,
        0x1289554b02b22083u64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0xb9feffffffffa9bbu64,
        0x1eabfffeb153ffffu64,
        0x6730d2a0f6b0f624u64,
        0x64774b84f38512bfu64,
        0x4b1ba7b6434bacd7u64,
        0x1a0111ea397fe69au64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x00000000000000f0u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&Fq2::zero());
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);

    let c0 = Fq::from_repr(FqRepr([
        0x0000000000076980u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x00000000003b4c00u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0xe24baa0a898b47e0u64,
        0x92afb1b88e09c84cu64,
        0xf16d677192b7b78au64,
        0xab1dd12189c47c0eu64,
        0xc30f74ce786d38e9u64,
        0x0cc49de633f05c98u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x936dda4aedcab1e1u64,
        0x08261a18f1038bdbu64,
        0x0c08dea79dde085du64,
        0x9002d76a3ed1ffd2u64,
        0x185ab763985ff885u64,
        0x00bab7cc25639665u64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x00000000000002d0u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
    ]))
    .unwrap();
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

    let m1 = {
        let mut tmp = Fq2::one();
        tmp.negate();
        tmp
    };
    let c0 = Fq::from_repr(FqRepr([
        0x0000000000076980u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x00000000003b4c00u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0xd7b355f5767462cbu64,
        0x8bfc4e46234a37b2u64,
        0x75c36b2f63f93e99u64,
        0xb9597a6369c096b0u64,
        0x880c32e7cade73edu64,
        0x0d3c7404058f8a01u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x269125b51234f8cau64,
        0x1685e5e5c0507424u64,
        0x5b27f3f958d2edc7u64,
        0xd474741ab4b312edu64,
        0x32c0f052aaebb451u64,
        0x19465a1e141c5035u64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x00000000000002d0u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
        0x0000000000000000u64,
    ]))
    .unwrap();
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
    let p = G2::osswu_map(&m1);
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);

    let u = Fq2 {
        c0: Fq::from_repr(FqRepr([
            0xab85e284be119b06u64,
            0xd66a22e1357a0a78u64,
            0xb2b1eebb23f4f1beu64,
            0x86b44d0ab7ba6c3cu64,
            0x074fb6220dae0f91u64,
            0x15a020c28c99d05cu64,
        ]))
        .unwrap(),
        c1: Fq::from_repr(FqRepr([
            0x7b66fff074efe46eu64,
            0xa06efb36880a24b7u64,
            0x7e29eca1f704feafu64,
            0xe059e38b408dd4ceu64,
            0x85d3318e078dfebau64,
            0x198bfdcafe694646u64,
        ]))
        .unwrap(),
    };
    let c0 = Fq::from_repr(FqRepr([
        0xd4076d7c0779c5e9u64,
        0x69512cb25235be96u64,
        0xd4528ddb1ae277cau64,
        0x9178b403ff1b0d02u64,
        0x2a9cf07eafb075c8u64,
        0x0cc5413d19d7ee31u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x1723c70ff5d8e43eu64,
        0x950084d6c95e0298u64,
        0x4c4e713218054a18u64,
        0x706ec7cb756425e2u64,
        0xa1688a3096f0ce29u64,
        0x1244bfa13bb40cddu64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x2b9e65c5a45f6284u64,
        0x562dfe86dc588870u64,
        0x8a69c4278050e9acu64,
        0x434de6b4c66f904bu64,
        0x421528f4e98f1a9fu64,
        0x0497daba05b79191u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xf69322eda18d7fabu64,
        0x1abf30aa6936363cu64,
        0x6c85ce4cf4e4ad4bu64,
        0x7fe4e86d8ed16fa9u64,
        0xd7c622658b988b3au64,
        0x054fd4cc6a4e6375u64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x0aa387c68cf7be7cu64,
        0xc03089a5cc09a4e6u64,
        0xe5edbccff2482a38u64,
        0xe3d6eb25462de571u64,
        0x4e335fdb1f778f4cu64,
        0x0971c6b317c83ef8u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xa6f4d7ba451b0f2fu64,
        0x72ef6cd1ab5147ecu64,
        0x41d4017d111ec2f5u64,
        0xe2e1ff096b87f7a0u64,
        0xa03be00ec98045c5u64,
        0x0794934ff7688600u64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&u);
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);

    let u = Fq2 {
        c0: Fq::from_repr(FqRepr([
            0xb1bd962564b60ec9u64,
            0x6356bf7e70c57bb2u64,
            0xc71fb10d9c9a5fa5u64,
            0x5c008d651ae8136eu64,
            0xfddca0c49bee1c3du64,
            0x037d15e2364d56f6u64,
        ]))
        .unwrap(),
        c1: Fq::from_repr(FqRepr([
            0xaccaa54b876c3ce5u64,
            0x391251fcd979c2dbu64,
            0x97b8b673a62e9cd5u64,
            0xd51248d4a164d299u64,
            0x6452efcd734f861eu64,
            0x0a73a64e9a0e3483u64,
        ]))
        .unwrap(),
    };
    let c0 = Fq::from_repr(FqRepr([
        0xc4ce89448e0f4fffu64,
        0x9a0035a789cfef4du64,
        0x4f48e1ff0d7132aeu64,
        0x07e84e66e9de7ef0u64,
        0x2102e039cb42c6b8u64,
        0x06fecf2041672340u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xe69438ef448ad850u64,
        0x954c0d297661fc82u64,
        0xbf358e867c59336fu64,
        0xfac3601dcbf862eau64,
        0xd9aceaf887525cdbu64,
        0x0e51b74f813dfd6au64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x59cd709b1b3b8072u64,
        0xe7443d98c2842710u64,
        0x534e94af70924cc7u64,
        0x60145fe41820a6dcu64,
        0x4a55bb8953ef1f50u64,
        0x0a2f5fead4ee6a5au64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xa265d8fea6249f2au64,
        0xa272322d41e8e6c2u64,
        0xefdf5789af18c2f7u64,
        0x5481e012dc9ab6c2u64,
        0xd5109fb7200322abu64,
        0x023a35721ac1dbe8u64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0xfa45c82ffbf31d6au64,
        0xfff4cd0628e44310u64,
        0x8065e297b5fca25fu64,
        0xe0ecec0fca950701u64,
        0x172acf37d5768b2fu64,
        0x156e2906d93e57b2u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xb8f2155d6d22d679u64,
        0xdcebf94b00ad9a1fu64,
        0x7af62548d12dcbfau64,
        0xbfd083985de108aau64,
        0xa9f2c47b9121ae28u64,
        0x0fc0cfaf7e05cae4u64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&u);
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);

    let u = Fq2 {
        c0: Fq::from_repr(FqRepr([
            0x3c3f2e8bf79d6585u64,
            0x4ab24c9d117c5d0au64,
            0x72dee8ec987fc836u64,
            0xc2d2114b8268d659u64,
            0xa2c21bf1bff581bdu64,
            0x0cd12e27b03be553u64,
        ]))
        .unwrap(),
        c1: Fq::from_repr(FqRepr([
            0xf3d78b949ecf9984u64,
            0xb782934b8b8e9a5bu64,
            0x038e927090ae25ecu64,
            0xbb732a7b2a94725bu64,
            0x5963aa64cdf7cf76u64,
            0x0476e812529bea48u64,
        ]))
        .unwrap(),
    };
    let c0 = Fq::from_repr(FqRepr([
        0xe3b587b17b1189eeu64,
        0x5dcd7a9fa7d313bcu64,
        0x4a69ab196b196b08u64,
        0xaf7800d84621add2u64,
        0xe23e0d36f4e250e6u64,
        0x03da0f83f05e669bu64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x2a4d0ce4472b5853u64,
        0x681b53be33dfe404u64,
        0x0726bf58ac3ac4e7u64,
        0xcfd610a249430aa8u64,
        0x92b6b223335d7cbeu64,
        0x0ab192613fc39903u64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0xee225d0910042c8du64,
        0xf2300b82cfc678d6u64,
        0xd9e86a46843133f4u64,
        0x55096f7659be0913u64,
        0xf45f04fe8542b01cu64,
        0x12df415bff83ceadu64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x07617e62c58173d3u64,
        0x8abb5a736f18e25eu64,
        0x710536695390588du64,
        0x074b734430b0fca7u64,
        0xd7e7f015840b6c0eu64,
        0x02671c2426b1ab44u64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x9f6b2cf30df6aec1u64,
        0x890ad07b8473fd9cu64,
        0x8cf00156bc32fc39u64,
        0x3878473de1ad29b6u64,
        0x1afd0d4ee306d7a6u64,
        0x06844c10fc7fd08au64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xbc91c4acf76320fbu64,
        0x2e6116dfcb361845u64,
        0x247c90216d613d76u64,
        0x4a1e607efc18c456u64,
        0x5237c926a44662d3u64,
        0x133fab21ff4f228fu64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&u);
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);

    let u = Fq2 {
        c0: Fq::from_repr(FqRepr([
            0x24ca8a83e995a3eau64,
            0x872b72d96f8346a8u64,
            0x7821255f7b61bea4u64,
            0x55cec1fd89d1e8d1u64,
            0x33114ee21773f115u64,
            0x01cd4733f8d1813cu64,
        ]))
        .unwrap(),
        c1: Fq::from_repr(FqRepr([
            0xc2a59c091a5fd6e0u64,
            0xda2b2811c76dc106u64,
            0xa6fb05a057525206u64,
            0x9715325bc1ff013au64,
            0x006ec8a9da0659c8u64,
            0x088855bf8353b4c7u64,
        ]))
        .unwrap(),
    };
    let c0 = Fq::from_repr(FqRepr([
        0x159fbfd21f8014f8u64,
        0xde6623e7dad77114u64,
        0x64d66f8f65b8ce6fu64,
        0x9fa138196080ad3eu64,
        0xaa3ff37e7c6626b9u64,
        0x03e68e6ff72199a8u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xec8440c8f359b118u64,
        0xc3b8a1e3b95f49b2u64,
        0x5fbc945a4c2f4cffu64,
        0xb02af47182709e7eu64,
        0x090a35d95de5079fu64,
        0x10d838cb45cbff85u64,
    ]))
    .unwrap();
    let xo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x2f602386c64c5820u64,
        0x0a453b197e9ad9a7u64,
        0xe1af19735c34a642u64,
        0xfb8943c756999f3fu64,
        0xa9bc72de80dcdc15u64,
        0x119a58b809dad83fu64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0x413ac6090656fd6au64,
        0xbf9a6d121c92b812u64,
        0x8fac010cb96c3dacu64,
        0xef7318cc80616a9du64,
        0x9d5de9c237dd356fu64,
        0x0d732bfd0e8c702bu64,
    ]))
    .unwrap();
    let yo = Fq2 { c0, c1 };
    let c0 = Fq::from_repr(FqRepr([
        0x62c5768f4975c6d6u64,
        0xbc5fb7b66dba95e6u64,
        0x539bcb35920ec1c9u64,
        0x4d90263bef559ad1u64,
        0x2464bf59a471c752u64,
        0x009207a5c9d0d734u64,
    ]))
    .unwrap();
    let c1 = Fq::from_repr(FqRepr([
        0xb916632d76be2983u64,
        0x4b01c074b64687fcu64,
        0x4e46bfdc0809abc7u64,
        0x86d345d9e596f5a1u64,
        0x03e64d0fcc4e7edbu64,
        0x1802bf2004793570u64,
    ]))
    .unwrap();
    let zo = Fq2 { c0, c1 };
    let p = G2::osswu_map(&u);
    let G2 { x, y, z } = &p;
    assert_eq!(x, &xo);
    assert_eq!(y, &yo);
    assert_eq!(z, &zo);

    let mut rng = thread_rng();
    for _ in 0..32 {
        let input = Fq2::rand(&mut rng);
        let p = G2::osswu_map(&input);
        let G2 { x, y, z } = &p;
        check_g2_prime(x, y, z);
    }
}
