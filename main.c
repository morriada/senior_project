#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <rtl-sdr.h>

static rtlsdr_dev_t *dev = NULL;

void usage(void)
{
  printf(
	 "TODO:Usage Info\n"
	 );
  exit(1);
}

int main(int argc, char** argv)
{
  int i, r, opt;
  uint32_t dev_index = 0;
  uint32_t i2c_addr = 0;
  uint8_t bias_on = 0;
  uint8_t i2c_value = 0;
  uint32_t sample_rate = 2400000;
  uint32_t center_freq = 151500000;
  uint32_t if_freq = 0;
  
  while((opt = getopt(argc, argv, "d:b:a:c:h?")) != -1) {
    switch(opt) {
    case 'd':
      dev_index = atoi(opt);
      break;
    case 'b':
      bias_on = atoi(opt);
      break;
    case 'a':
      i2c_addr = strtol(opt, NULL, 16);
      printf("i2c address: %02x \n", i2c_addr);
      break;
    case 'v':
      i2c_value = strtol(opt, NULL, 16); 
      break;
    default:
      usage();
      break;
    }
  }
  
  r = rtlsdr_open(&dev, dev_index);

  // Set the sample rate of the rtl-sdr
  CHECK1(rtlsdr_set_sample_rate(dev, sample_rate));
  // Set dithering 0
  CHECK1(rtlsdr_set_dithering(dev, 0));
  // Set the IF frequency
  CHECK1(rtlsdr_set_if_freq(dev, if_freq));
  // Set the center frequency
  CHECK1(rtlsdr_set_center_freq(dev, center_freq));
  // Set the tuner gain mode to automatic
  CHECK1(rtlsdr_set_tuner_gain_mode(dev, 0));
  // Set the bias tee by setting the gpio bit 0 to bias_off
  CHECK1(rtlsdr_set_bias_tee(dev, bias_noise));
  // Set rtlsdr repeater for the i2communication via RTL2838
  CHECK1(rtlsdr_set_i2c_repeater(dev, 1));
  // Set register to the output
  CHECK1(rtlsdr_i2c_write_reg(dev, i2c_addr, 0x03, 00));
  // Set value to the register as described in the table
  CHECK1(rtlsdr_i2c_write_reg(dev, i2c_addr, 0x01, i2c_value));
  // Close the i2c_repeater
  CHECK1(rtlsdr_set_i2c_repeater(dev, 0));
  // Reset the buffer
  CHECK1(rtlsdr_reset_buffer(dev));

  return 0;
}
