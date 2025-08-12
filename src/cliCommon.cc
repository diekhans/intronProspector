/* cliCommon.cc - common CLI parsing functions.

    Copyright (c) 2018-2025, Mark Diekhans, University of California, Santa Cruz

This code is derived from the regtools package available at:

  https://github.com/griffithlab/regtools

We thank them for providing this software under a flexible license:

    Copyright (c) 2015, The Griffith Lab

    Author: Avinash Ramu <aramu@genome.wustl.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include "cliCommon.hh"
#include "version.hh"


static struct option common_long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {"junction-bed", required_argument, NULL, 'j'},
    {"intron-bed", required_argument, NULL, 'n'},
    {"intron-bed6", required_argument, NULL, 'b'},
    {"intron-calls", required_argument, NULL, 'c'},
    {NULL, 0, NULL, 0}
};

static const char *common_short_options = "hvj:n:b:c:";

/* count up to NULL */
static unsigned count_long_options(struct option* long_options) {
    unsigned cnt = 0;
    for (; long_options[cnt].name != NULL; cnt++) {
        continue;
    }
    return cnt;
}

/* merge common long options with program-specific */
struct option *add_common_long_options(struct option* long_options) {
    static struct option end_marker = {NULL, 0, NULL, 0};
    unsigned num_opts = count_long_options(common_long_options) +
        count_long_options(long_options);

    struct option* combined_options = new struct option[num_opts + 1];
    int iout = 0;
    for (int i = 0; common_long_options[i].name != NULL; i++) {
        combined_options[iout++] = common_long_options[i];
    }
    for (int i = 0; long_options[i].name != NULL; i++) {
        combined_options[iout++] = long_options[i];
    }
    combined_options[iout] = end_marker;
    return combined_options;
}


/* merge common short options with program-specific */
string add_common_short_options(const string& long_options) {
    return  common_short_options + long_options;
}
