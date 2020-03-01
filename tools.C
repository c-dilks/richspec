// convert maroc channel to pmt pixel number
Int_t chan2pix(Int_t chan) {
  Int_t c2p[] = {
    60, 58, 59, 57, 52, 50, 51, 49,
    44, 42, 43, 41, 36, 34, 35, 33,
    28, 26, 27, 25, 20, 18, 19, 17,
    12, 10, 11,  9,  4,  2, 3,   1,
     5,  7,  6,  8, 13, 15, 14, 16,
    21, 23, 22, 24, 29, 31, 30, 32,
    37, 39, 38, 40, 45, 47, 46, 48,
    53, 55, 54, 56, 61, 63, 62, 64};
  return c2p[chan%64];
};

// convert pixel number to row and column
Int_t xPix(Int_t pix) { return (pix-1) % 8; };
Int_t yPix(Int_t pix) { return 7 - (pix-1) / 8; };

// convert maroc channel to pmt number
Int_t chan2pmt(Int_t chan) { return chan/64; };


