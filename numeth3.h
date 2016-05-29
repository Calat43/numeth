#pragma once

#include <vector>
#include "ploter.h"
#include <iostream>

void maccormack( AnimPloter & ro_ploter_win,
                 AnimPloter & u_ploter_win,
                 AnimPloter & p_ploter_win,
                 std::ostream & fout );
