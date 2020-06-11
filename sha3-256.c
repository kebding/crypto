/* sha3-256.c

   This program reads an input file and produces its SHA3-256 hash.
   The hash produced is printed to the terminal.
   This program assumes the input to be byte-aligned, i.e. the number of bits
   in the file is divisible by 8.

   Author: Kyle Ebding
   Created on 2020-04-16
*/

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#define KECCAK_ROUNDS 24
#define OUTPUT_BITS 256
#define OUTPUT_BYTES (OUTPUT_BITS / 8)  // 32
#define RATE 1088
#define RATE_BYTES (RATE / 8)    // 136
#define CAPACITY 512
#define CAPACITY_BYTES (CAPACITY / 8)  // 64

#define ROTL64(x, r) (((x) << (r)) | ((x) >> (64 - (r))))
#define ROTR64(x, r) (((x) >> (r)) | ((x) << (64 - (r))))

void keccak(uint64_t *state[5], int verbose);
void pad_input(uint8_t* input, int input_len);
void absorb_input(uint8_t *input, uint64_t *state[5]);
void print_state(uint64_t *state[5]);
void print_input(uint8_t *input);
void print_lanes(uint64_t *lanes[5]);
void print_3d(uint64_t *state[5]);

int main(int argc, char** argv) {

    if (argc < 2) {
        printf("USAGE: ./sha3 <file_to_hash>\n");
        return -1;
    }

    // open the file to hash
    FILE *infile;
    infile = fopen(argv[1], "rb");
    if (infile == NULL) {
        printf("ERROR: could not read file %s\n", argv[1]);
        return -1;
    }

    int verbose = 0;
    if (argc > 2 && argv[2][0] == '-' && argv[2][1] == 'v') {
        verbose = 1;
    }

    int i, j, bytes_read;
    uint8_t *input = (uint8_t*)calloc(RATE_BYTES + CAPACITY_BYTES, sizeof(uint8_t));
    uint64_t *state[5];
    for (i = 0; i < 5; i++) {
        state[i] = (uint64_t*)calloc(5, sizeof(uint64_t));
    }

    bytes_read = fread(input, 1, RATE_BYTES, infile);
    while(1) {
        if (verbose) {
            printf("About to absorb data\n");
            printf("State (in bytes)\n");
            print_state(state);
        }
        pad_input(input, bytes_read);
        if (verbose) {
            printf("Data to be absorbed\n");
            print_input(input);
        }
        absorb_input(input, state);
        if (verbose) {
            printf("XOR'd state (in bytes)\n");
            print_state(state);
            printf("XOR'd state (as lanes of integers)\n");
            print_lanes(state);
        }
        keccak(state, verbose);
        if (bytes_read == RATE_BYTES) { // i.e. if there is more input to read
            memset(input, 0, RATE_BYTES + CAPACITY_BYTES);
            bytes_read = fread(input, 1, RATE_BYTES, infile);
        }
        else break;
    }

    if (verbose) {
        printf("After permutation\n");
        print_state(state);
        printf("state (as lanes of integers)\n");
        print_lanes(state);
        printf("Hash val is\n");
        // prints first 256 bits of state as bytes in two rows
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 8; j++) {
                printf("%02X ", (uint8_t)(state[0][i] >> (j * 8)));
            }
            if ((i + 1) % 2 == 0) {
                printf("\n");
            }
        }
    }
    // print the result
    printf("0x%lx%lx%lx%lx\n", state[0][0], state[0][1], state[0][2],
            state[0][4]);

    for (i = 0; i < 5; i++)
        free(state[i]);
    free(input);
    // close the input file
    fclose(infile);
}


/* pad_input

   This function pads the input according to the SHA3-256 spec.
   Note that the parameter input is known to be 200 bytes and that the
   parameter input_len is to know where to begin padding.
*/
void pad_input(uint8_t* input, int input_len) {
    int i, q;
    if (input_len < RATE_BYTES) {
        // pad the input
        q = RATE_BYTES - (input_len % RATE_BYTES);
        if (q == 1) {   // pad with 0x86
            input[input_len] = 0x86;
        }
        else if (q == 2) {  // pad with 0x0680
            input[input_len] = 0x06;
            input[input_len + 1] = 0x80;
        }
        else {  // pad with 0x06 [q-2 zero bytes] 0x80
            input[input_len] = 0x06;
            for (i = 1; i < q - 2; i++) {
                input[input_len + i] = 0x00;
            }
            input[input_len + i + 1] = 0x80;
        }
    }
}

/* absorb_input

   This function takes a block of the input message and absorbs it into the
   state in preparation for applying Keccak to the state.
*/
void absorb_input(uint8_t *input, uint64_t *state[5]) {

    int i, j;

    // XOR the input into the state lane-by-lane
    uint64_t input_word;
    for (i = 0; i < RATE_BYTES / 8; i++) {
        input_word = 0;
        for (j = 7; j >= 0; j--) {
            input_word <<= 8;
            input_word |= input[(i * 8) + j];
        }
        state[(i / 5)][(i % 5)] ^= input_word;
    }

}


/* keccak

   This function performs Keccak-f[1600] on the state.
   This is based on the pseudo-code provided at
        https://keccak.team/keccak_specs_summary.html
   The rho/pi step is from
        https://github.com/mjosaarinen/tiny_sha3/blob/master/sha3.c
*/
void keccak(uint64_t *state[5], int verbose) {

    int i, j, round;
    // need to use heap memory to prevent stack buffer overflow
    uint64_t temp_lane;
    uint64_t *temp_row;
    temp_row = (uint64_t*)calloc(5, sizeof(uint64_t));

    const uint64_t round_constants[24] = {
        0x0000000000000001, 0x0000000000008082, 0x800000000000808a,
        0x8000000080008000, 0x000000000000808b, 0x0000000080000001,
        0x8000000080008081, 0x8000000000008009, 0x000000000000008a,
        0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
        0x000000008000808b, 0x800000000000008b, 0x8000000000008089,
        0x8000000000008003, 0x8000000000008002, 0x8000000000000080,
        0x000000000000800a, 0x800000008000000a, 0x8000000080008081,
        0x8000000000008080, 0x0000000080000001, 0x8000000080008008
    };
    const int rho_rotations[24] = {
         1,  3,  6, 10, 15, 21, 28, 36, 45, 55,  2, 14,
        27, 41, 56,  8, 25, 43, 62, 18, 39, 61, 20, 44
    };
    const int pi_shifts[24] = {
        10, 7,  11, 17, 18, 3, 5,  16, 8,  21, 24, 4,
        15, 23, 19, 13, 12, 2, 20, 14, 22, 9,  6,  1
    };

    for (round = 0; round < KECCAK_ROUNDS; round++) {

        if (verbose) {
            printf("Round #%d\n", round);
            //print_3d(state);
        }
        // theta
        // get the parity of each column
        for (i = 0; i < 5; i++) {
            temp_row[i] = 0;
            for (j = 0; j < 5; j++) {
                temp_row[i] ^= state[j][i];
            }
        }
        // get the XOR sum of the temp_row of neighboring columns
        for (i = 0; i < 5; i++) {
            // index (i+4)%5 is i-1 and index (i+1)%5 is i+1 (circularly)
            temp_lane = temp_row[(i + 4) % 5] ^ ROTL64(temp_row[(i + 1) % 5], 1);
            // add column's neighbors' temp_row' sum into the column
            for (j = 0; j < 5; j++) {
                state[j][i] ^= temp_lane;
            }
        }
        if (verbose) {
            printf("After Theta\n");
            print_state(state);
            //print_3d(state);
        }

        // rho and pi
        temp_lane = state[0][1];
        for (i = 0; i < 24; i++) {
            j = pi_shifts[i];
            temp_row[0] = state[j / 5][j % 5];
            state[j / 5][j % 5] = ROTL64(temp_lane, rho_rotations[i]);
            temp_lane = temp_row[0];
        }

        if (verbose) {
            printf("After Rho and Pi\n");
            print_state(state);
            //print_3d(state);
        }

        // chi
        for (i = 0; i < 5; i++) {
            // copy the row
            for (j = 0; j < 5; j++) {
                temp_row[j] = state[i][j];
            }
            // if the next two lanes to the right (circularly) exhibit the
            // pattern '01' then flip the current lane (for each bit in the
            // lane)
            for (j = 0; j < 5; j++) {
                state[i][j] ^=
                    ((~temp_row[(j + 1) % 5]) & temp_row[(j + 2) % 5]);
            }
        }
        if (verbose) {
            printf("After Chi\n");
            print_state(state);
            //print_3d(state);
        }

        // iota
        // break symmetry by mixing in the round constant
        state[0][0] ^= round_constants[round];
        if (verbose) {
            printf("After Iota\n");
            print_state(state);
            //print_3d(state);
        }
    }
    free(temp_row);
    return;
}


/* print_state

   This function prints out the state data byte-by-byte.
   This will be used when the script is in verbose mode to print the
   intermediate state of the data for testing. The output format is designed
   to match the format of the official test vectors (see README).

   This function assumes that the input in a 5x5 array of long unsigned ints.
*/
void print_state(uint64_t *state[5]) {
    int i, j, k;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++) {
            for (k = 0; k < 8; k++) {
                printf("%02X ", (uint8_t)(state[i][j] >> (k * 8)));
            }
            if ((i + j) % 2 == 1) {
                printf("\n");
            }
        }
    }
    printf("\n");
}

/* print_input

   This function prints the input for the round of SHA3 byte-by-byte.
   This will be used when the script is in verbose mode.
   The output format is designed to match the format of the official test
   vectors (see README).

   This function assumes that the input is a byte array of length RATE_BYTES +
   CAPACITY_BYTES
*/
void print_input(uint8_t *input) {
    int i;
    for (i = 0; i < RATE_BYTES + CAPACITY_BYTES; i++) {
        printf("%02X ", input[i]);
        if ((i + 1) % 16 == 0) {
            printf("\n");
        }
    }
    printf("\n");
}

/* print_lanes

   This function prints the lanes of data as integers, as seen in the official
   test vectors. This function will be used in verbose mode to print
   intermediate values for testing.
   The rows/columns are reversed in the print vs the data storage because the
   official test vectors use column-major indexing.
*/
void print_lanes(uint64_t *lanes[5]) {
    int i, j;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++) {
            printf("  [%d, %d] = %016lx\n", j, i, lanes[i][j]);
        }
    }
}

/* print_3d

   This function prints the lanes in row/column order. While not truly 3D, it
   does show each byte in a way that can easily be examined alongside the
   byte's neighbors in the structure for easier bit comparison than in the
   format given by print_lanes or print_state.
   This was used for debugging only and is now unused.
*/
void print_3d(uint64_t *state[5]) {
    int i;
    for (i = 0; i < 5; i++) {
        printf("row %d: %016lx %016lx %016lx %016lx %016lx\n", i,
                state[i][0], state[i][1], state[i][2], state[i][3], state[i][4]);
    }
}
