//
// A simple example to get you started with the library.
// You can compile and run this example like so:
//
//   make example
//   ./example
//
//  Warning: If your compiler does not fully support C++11, some of
//  this example may require changes.
//

#include "headers/codecfactory.h"
#include "headers/deltautil.h"
#include "headers/synthetic.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define BATCH_SIZE 8
#define LG_BATCH_SIZE 3
#define LG_BATCH_SIZE_PLUS_ONE 4
#define WORD_SIZE 32
#define LG_WORD_SIZE 5
#define NUM_WIDTHS 16
#define LG_NUM_WIDTHS 4

typedef uint8_t uint1_t;

typedef struct {
  uint1_t fixed_cheaper;
  uint16_t cost;
} cost_info;

static inline uint8_t bit_length(uint32_t word) {
  return word == 0 ? 1 : WORD_SIZE - __builtin_clz(word);
}

static inline uint32_t bit_select(uint32_t word, uint32_t upper, uint32_t lower) {
  uint32_t num_bits = upper - lower + 1;
  return (word >> lower) & ((1L << num_bits) - 1);
}

static inline cost_info compute_cost(uint8_t width, uint16_t bit_count[3]) {
  uint8_t num_exceptions = BATCH_SIZE - bit_count[0];
  uint8_t fixed_cost = LG_WORD_SIZE + num_exceptions * bit_count[1];
  uint8_t varint_cost = bit_count[2]; // will be 0 if there are no exceptions
  uint8_t common_exception_cost = (num_exceptions > 0 ? 1 : 0) + num_exceptions * MAX(LG_BATCH_SIZE, 1);
  return (cost_info) {fixed_cost <= varint_cost, (uint16_t)(width * bit_count[0] +
    common_exception_cost + (fixed_cost <= varint_cost ? fixed_cost : varint_cost))};
}

void run(uint32_t *input, uint32_t input_count, uint32_t *output, uint32_t *output_count) {
  uint32_t input_idx = 0;
  uint32_t output_idx = 0;
  uint32_t out_buf = 0;
  uint8_t out_buf_bits = 0;
  uint8_t bits_to_varint_len[33];
  uint8_t bit_widths[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 20, 32};
  uint16_t bit_counts[NUM_WIDTHS][3];
  uint8_t words_consumed = 0;
  uint32_t buffer[BATCH_SIZE];

  for (int i = 1; i < 33; i++) {
    bits_to_varint_len[i] = (i + 7 - 1) / 7 * 8;
  }
  for (int i = 0; i < NUM_WIDTHS; i++) {
    bit_counts[i][0] = 0;
    bit_counts[i][1] = 0;
    bit_counts[i][2] = 0;
  }

  for (uint32_t i = input_idx; i < input_count; i++) {
    #define ADD_TO_OUTPUT(w0, b0) do {\
              uint32_t word = w0;\
              uint8_t num_bits = b0;\
              if (num_bits + out_buf_bits >= 32) {\
                output[output_idx++] = ((word & ((1L << (32 - out_buf_bits)) - 1)) << out_buf_bits) | out_buf;\
                out_buf = word >> (32 - out_buf_bits);\
                out_buf_bits = num_bits + out_buf_bits - 32;\
              } else {\
                out_buf = (word << out_buf_bits) | out_buf;\
                out_buf_bits += num_bits;\
              }\
            } while (0)
    // could replace this with the direct computatation instead of table lookup
    #define VARINT_LEN(b0) bits_to_varint_len[b0]

    uint8_t cur_bit_length = bit_length(input[i]);
    for (int j = 0; j < NUM_WIDTHS; j++) {
      if (cur_bit_length <= bit_widths[j]) {
        bit_counts[j][0] += 1;
      } else {
        if (cur_bit_length > bit_counts[j][1]) {
          bit_counts[j][1] = cur_bit_length;
        }
        bit_counts[j][2] += VARINT_LEN(cur_bit_length);
      }
    }
    buffer[words_consumed++] = input[i];
    if (words_consumed == BATCH_SIZE) {
      uint8_t min_width_idx = 0;
      cost_info min_cost = {0, 65535}; // assumes that 65535 is an unreachable cost value
      for (int j = 0; j < NUM_WIDTHS; j++) {
        cost_info cur_cost = compute_cost(bit_widths[j], bit_counts[j]);
        if (cur_cost.cost < min_cost.cost) {
          min_width_idx = j;
	  min_cost = cur_cost;
        }
      }
      ADD_TO_OUTPUT(min_width_idx, MAX(LG_NUM_WIDTHS, 1));
      ADD_TO_OUTPUT(BATCH_SIZE - bit_counts[min_width_idx][0], LG_BATCH_SIZE_PLUS_ONE);
      for (int j = 0; j < BATCH_SIZE; j++) {
        if (bit_length(buffer[j]) <= bit_widths[min_width_idx]) {
          ADD_TO_OUTPUT(buffer[j], bit_widths[min_width_idx]);
        }
      }
      if (bit_counts[min_width_idx][0] < BATCH_SIZE) {
        if (min_cost.fixed_cheaper) {
          // fixed exceptions
          ADD_TO_OUTPUT(0, 1);
          ADD_TO_OUTPUT(bit_counts[min_width_idx][1] - 1, LG_WORD_SIZE);
        } else {
          // varint exceptions
          ADD_TO_OUTPUT(1, 1);
        }
        for (int j = 0; j < BATCH_SIZE; j++) {
          uint8_t word_bit_length = bit_length(buffer[j]);
          if (word_bit_length > bit_widths[min_width_idx]) {
            ADD_TO_OUTPUT(j, MAX(LG_BATCH_SIZE, 1));
            if (min_cost.fixed_cheaper) { // fixed
              ADD_TO_OUTPUT(buffer[j], bit_counts[min_width_idx][1]);
            } else {
              for (int k = 0; k < word_bit_length; k += 7) {
                if (k + 7 < word_bit_length) {
                  ADD_TO_OUTPUT(bit_select(buffer[j], k + 6, k) | (1 << 7), 8);
                } else {
                  ADD_TO_OUTPUT(bit_select(buffer[j], word_bit_length - 1, k), 8);
                }
              }
            }
          }
        }
      }
      for (int j = 0; j < NUM_WIDTHS; j++) {
        bit_counts[j][0] = 0;
        bit_counts[j][1] = 0;
        bit_counts[j][2] = 0;
      }
      words_consumed = 0;
    }
  }
  if (out_buf_bits > 0) {
    output[output_idx++] = out_buf;
  }
  *output_count = output_idx;
}

int main() {
  using namespace FastPForLib;

  // We pick a CODEC
  IntegerCODEC &codec = *CODECFactory::getFromName("optpfor");
  // could use others, e.g., "simdbinarypacking", "varintg8iu", "simdfastpfor256"
  ////////////
  //
  // create a container with some integers in it
  //
  // for this example, we will not assume that the
  // integers are in sorted order
  //
  // (Note: You don't need to use a vector.)
  //
  UniformDataGenerator clu;
  uint32_t N = 1U << 25;
  auto mydata = clu.generateUniform(N, 1U << 29);
  /*
  size_t N = 10 * 1000;
  std::vector<uint32_t> mydata(N);
  for (uint32_t i = 0; i < N; i += 150)
    mydata[i] = i;
  */
  Delta::deltaSIMD(mydata.data(), mydata.size());
  //
  // the vector mydata could contain anything, really
  //
  ///////////
  //
  // You need some "output" container. You are responsible
  // for allocating enough memory.
  //
  std::vector<uint32_t> compressed_output(N + 1024);
  // N+1024 should be plenty
  //
  //
  size_t compressedsize = compressed_output.size();
  codec.encodeArray(mydata.data(), mydata.size(), compressed_output.data(),
                    compressedsize);
  //
  // if desired, shrink back the array:
  compressed_output.resize(compressedsize);
  compressed_output.shrink_to_fit();
  // display compression rate:
  std::cout << std::setprecision(3);
  std::cout << "optpfor "
            << 32.0 * static_cast<double>(compressed_output.size()) /
                   static_cast<double>(mydata.size())
            << " bits per integer. " << std::endl;

  uint32_t output_count;
  compressed_output.resize(N + 1024);
  run(mydata.data(), (uint32_t)mydata.size(), compressed_output.data(), &output_count);
  std::cout << "our approach "
            << 32.0 * static_cast<double>(output_count) /
                   static_cast<double>(mydata.size())
            << " bits per integer. " << std::endl;
  /*
  //
  // You are done!... with the compression...
  //
  ///
  // decompressing is also easy:
  //
  std::vector<uint32_t> mydataback(N);
  size_t recoveredsize = mydataback.size();
  //
  codec.decodeArray(compressed_output.data(), compressed_output.size(),
                    mydataback.data(), recoveredsize);
  mydataback.resize(recoveredsize);
  //
  // That's it!
  //
  if (mydataback != mydata)
    throw std::runtime_error("bug!");

  // If you need to use differential coding, you can use
  // calls like these to get the deltas and recover the original
  // data from the deltas:
  Delta::deltaSIMD(mydata.data(), mydata.size());
  Delta::inverseDeltaSIMD(mydata.data(), mydata.size());
  // be mindful of CPU caching issues

  // If you do use differential coding a lot, you might want 
  // to check out these other libraries...
  // https://github.com/lemire/FastDifferentialCoding
  // and
  // https://github.com/lemire/SIMDCompressionAndIntersection
  */
}
