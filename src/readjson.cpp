/***********************************************************************************************
 * Copyright (C) 2017 Qinsi Wang.  All rights reserved.
 * By using this software the USER indicates that he or she has read, understood and will comply
 * with the following:
 *
 * 1. The USER is hereby granted non-exclusive permission to use, copy and/or
 * modify this software for internal, non-commercial, research purposes only. Any
 * distribution, including commercial sale or license, of this software, copies of
 * the software, its associated documentation and/or modifications of either is
 * strictly prohibited without the prior consent of the authors. Title to copyright
 * to this software and its associated documentation shall at all times remain with
 * the authors. Appropriated copyright notice shall be placed on all software
 * copies, and a complete copy of this notice shall be included in all copies of
 * the associated documentation. No right is granted to use in advertising,
 * publicity or otherwise any trademark, service mark, or the name of the authors.
 *
 * 2. This software and any associated documentation is provided "as is".
 *
 * THE AUTHORS MAKE NO REPRESENTATIONS OR WARRANTIES, EXPRESSED OR IMPLIED,
 * INCLUDING THOSE OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
 * USE OF THE SOFTWARE, MODIFICATIONS, OR ASSOCIATED DOCUMENTATION WILL NOT
 * INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER INTELLECTUAL PROPERTY
 * RIGHTS OF A THIRD PARTY.
 *
 * The authors shall not be liable under any circumstances for any direct,
 * indirect, special, incidental, or consequential damages with respect to any
 * claim by USER or any third party on account of or arising from the use, or
 * inability to use, this software or its associated documentation, even if the
 * authors have been advised of the possibility of those damages.
 * ***********************************************************************************************/


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "readjson.hpp"

using std::string;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::vector;

std::vector<int> readjson (string const & rulenameA, string const & rulenameB, string const & influenceType, int simulationround) {

    //int ret = 1;

    string line;
    string line2;
    ifstream fluxtable ("activities.json");
    vector<string> ruleset;
    vector<string> impactset;
    //string ruleAnameinput = rulenameA;
    string ruleAname = "\"" + string(rulenameA) + "\"";
    //string ruleBnameinput = rulenameB;
    string ruleBname = "\"" + string(rulenameB) + "\"";
    string delim = ","; // for rules
    string rulesline;
    float accImpact = 0.0;
    int impactfinal;
    int ruleAocc = 0;
    vector<int> returnres;
    

    ofstream fluxfile ("fluxfile_" + std::to_string(simulationround) + ".out");
    //ofstream fluxcheking ("fluxchecking.txt");

    if (fluxtable.is_open())
    {
    
        while (getline (fluxtable, line))
        {

            if (line.substr(0, 1) == "r")
            {
                rulesline = line;
                break;
            }
        }
    
        string rulesnames = rulesline.substr(7, rulesline.length()-9);
    
        size_t prev = 0, pos = 0;
    
        do
        {
            pos = rulesnames.find(delim, prev);
            if (pos == string::npos) pos = rulesnames.length();
            string rule = rulesnames.substr(prev, pos-prev);
            if (!rule.empty()) {
                ruleset.push_back(rule);
            }
            prev = pos + delim.length();
        } while (pos < rulesnames.length() && prev < rulesnames.length());
      
        int ruleAind = std::find(ruleset.begin(), ruleset.end(), ruleAname) - ruleset.begin();
        int ruleBind = std::find(ruleset.begin(), ruleset.end(), ruleBname) - ruleset.begin();
        
        string ruleAindex = std::to_string(ruleAind);
        string ruleBindex = std::to_string(ruleBind);
        
        
        string dataline;
        string impact;
        //int linenum = 0;
        
        while (getline (fluxtable, line2))
        {
            prev = 0;
            pos = 0;
            //linenum += 1;
            
            if (line2.substr(0, 1) == "[")
            {
                dataline = line2.substr(1, line2.length()-3);
                
                do
                {
                    pos = dataline.find(delim, prev);
                    if (pos == string::npos) pos = dataline.length();
                    impact = dataline.substr(prev, pos-prev);
                    
                    if (!impact.empty()) {
                        impactset.push_back(impact);
                    }
                    prev = pos + delim.length();
                } while (pos < dataline.length() && prev < dataline.length());

                if (impactset[0] == ruleAindex)
                {
                    //ruleAocc = linenum;
                    ruleAocc += 1;
                    accImpact += std::stof(impactset[ruleBind]);
                    fluxfile << accImpact << '\n';
                    /*
                    if (accImpact < 0){
                        fluxcheking << "1" << '\n';
                    }else{
                        fluxcheking << "0" << '\n';
                    }
                    */
                }
        
                impactset.clear();
            }
        }
        
        fluxtable.close();
        fluxfile.close();
        //fluxcheking.close();
        //cout << ruleAocc << endl;

    }else {
      cout << "Unable to open the given flux file";
      exit (EXIT_FAILURE);
    }
    
    std::string s1 = "positive";
    std::string s2 = "negative";
    std::string influType = string(influenceType);
    
    if (influType == s1 && accImpact > 0) {
        impactfinal = 1;
    } else if (influType == s2 && accImpact < 0) {
        impactfinal = 1;
    } else {impactfinal = 0;}
    
    returnres.push_back(impactfinal);
    returnres.push_back(ruleAocc);
    return returnres;
}
