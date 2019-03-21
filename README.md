# Pairing-fork

TL;DR: This is a fork of the [pairing library](https://github.com/zkcrypto/pairing), with additional features that make the library easy to use for BLS signature scheme.
The exact GitDiff is available [here](https://github.com/algorand/pairing-fork/compare/183a64b08e9dc7067f78624ec161371f1829623e...master)


## Details of changes

### Dependencies

* SHA2 library

### Additional functions

* hash to E1/E2
   * Based on SHA384
   * Try-and-increment method

* hash to G1/E2
   * hash to E1/E2 than multiply by co-factor

* Montgomery multiplications for E1/E2
   * number of doublings and additions are constant regardless the scalar

* Shamir's trick: simultaneous 2 multiplications

### Wrappers of existing functions

* cast strings to and from E1/E2

* membership testing

* two pairing product operation
