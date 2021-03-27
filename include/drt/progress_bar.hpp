#ifndef DRT_PROGRESS_BAR_HPP
#define DRT_PROGRESS_BAR_HPP
#include <iostream>
#include <chrono>

namespace drt {

void display_progress(double progress)
{
    int barWidth = 100;
    std::cout << "\033[0;32m[";
    int pos = static_cast<int>(barWidth * progress);
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "%\033[0m\r";
    std::cout.flush();
}

}
#endif
