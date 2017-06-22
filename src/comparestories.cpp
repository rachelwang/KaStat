#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include "comparestories.hpp"

using namespace std;

int comparestories(vector<string> & storyA, vector<string> & storyB) {
    
    std::vector<string>::iterator position;
    
    if (storyA.size() != storyB.size()) {
        return 0;
    } else {
        for (unsigned storyAi = 0; storyAi < storyA.size(); storyAi++) {
            position = std::find(storyB.begin(), storyB.end(), storyA[storyAi]);
            if (position != storyB.end()) {
                storyB.erase(position);
            } else {
                return 0;
                break;
            }
        }
        return 1;
    }
}
