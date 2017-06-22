#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include "readdot.hpp"

using namespace std;

std::vector<string> readdot(string const & storyfile) {
    
    string line;
    ifstream dotfile;
    string filename = storyfile;
    // ofstream edgefile;
    vector<string> edgefile;
    vector<string> eventind;
    vector<string> eventnam;
    vector<string> edgeset;
    vector<string> sourceset;
    vector<string> targetset;
    vector<string> reducededgeset;
    string reducededge;
    
    dotfile.open (filename.c_str());
    if (dotfile.fail()) {
        cout << "\nInvaild file name.\n";
        exit (EXIT_FAILURE);
    }
    //edgefile.open(("out_"+filename+".txt").c_str());
    
    std::ifstream readFile(filename);
    
    size_t pos = 0;
    while (getline (readFile, line))
    {
        if (line.find("[label=\"") != std::string::npos)
        {
            std::string start_delim = "[label=\"";
            std::string stop_delim = "\",";
            unsigned first_delim_pos = line.find(start_delim);
            unsigned end_pos_of_first_delim = first_delim_pos + start_delim.length();
            unsigned last_delim_pos = line.find(stop_delim);
            std::string eventname = line.substr(end_pos_of_first_delim,
                            last_delim_pos - end_pos_of_first_delim);
            line.erase(std::remove(line.begin(), line.end(), ' '),
                       line.end());
            pos = line.find("[");
            eventind.push_back(line.substr(0, pos));
            eventnam.push_back(eventname);
            //edgefile << line.substr(0, pos) << " " << eventname << "\n";
        }
    }
    
    readFile.clear();
    readFile.seekg(0, ios::beg);
    string source, target;
    string nusource, nutarget;
    
    while (getline (readFile, line))
    {
        if (line.find("->") != std::string::npos)
        {
            while(line.find("[") != std::string::npos) {
                size_t Beg = line.find("[");
                line.erase(Beg, line.find("\n", Beg) - Beg);
            }
            line.erase(std::remove(line.begin(), line.end(), ' '),
                       line.end());
            if (std::find(edgeset.begin(), edgeset.end(), line) == edgeset.end()) {
                edgeset.push_back(line);
                //edgefile << line << "\n";
                
                pos = line.find("->");
//                source = line.substr(0, pos);
//                target = line.substr(pos+2, line.length()-pos-2);
//                for ( unsigned j = 0; j < eventind.size(); j++) {
//                    if (eventind[j] == source) nusource = eventnam[j];
//                    if (eventind[j] == target) nutarget = eventnam[j];
//                }
                //edgefile << nusource << " " << nutarget << "\n";
                //edgefile << line.substr(0, pos) << " " << line.substr(pos+2, line.length()-pos-2) << "\n";
                
                sourceset.push_back(line.substr(0, pos));
                targetset.push_back(line.substr(pos+2, line.length()-pos-2));
                
                
            }
        }
    }
    
    // for each node, find all reachable nodes from it
    vector<string> reachset4node;
    vector<string> reachset;
    vector<string> allreachset;
    string reachnodes;
    string currnode;
    
    for (unsigned nodei = 0; nodei < eventind.size(); nodei++) {
        
        
        if (std::find(reachset.begin(), reachset.end(), eventind[nodei]) == reachset.end()) {
            reachset4node.push_back(eventind[nodei]);
            reachset.push_back(eventind[nodei]);
        }
        reachnodes = eventind[nodei];
        while (!reachset4node.empty()) {
            currnode = reachset4node[0];
            reachset4node.erase(reachset4node.begin());
            for (unsigned edgei = 0; edgei < sourceset.size(); edgei++) {
                if (sourceset[edgei] == currnode && std::find(reachset.begin(), reachset.end(), targetset[edgei]) == reachset.end()) {
                    reachset4node.push_back(targetset[edgei]);
                    reachset.push_back(targetset[edgei]);
                    reachnodes = reachnodes + " " + targetset[edgei] + " ";
                }
            }
        }
        allreachset.push_back(reachnodes);
        //cout << reachnodes << endl;
        reachset.clear();
    }
    
    // given two nodes, if their names are the same, and there is no path between them, we can say that they are equal
    for (unsigned eventi = 0; eventi < eventind.size(); eventi++) {
        for (unsigned eventj = eventi; eventj < eventind.size(); eventj++) {
            
            if (eventnam[eventi] == eventnam[eventj] && allreachset[eventi].find(" " + eventind[eventj] + " ") == std::string::npos && allreachset[eventj].find(" " + eventind[eventi] + " ") == std::string::npos && eventind[eventi] != eventind[eventj]) {
                //cout << eventind[eventj] << " " << eventind[eventi] << endl;
                std::replace (sourceset.begin(), sourceset.end(), eventind[eventj], eventind[eventi]);
                std::replace (targetset.begin(), targetset.end(), eventind[eventj], eventind[eventi]);
            }
        }
    }
    
    // remove duplicated edges
    for (unsigned edgej = 0; edgej < sourceset.size(); edgej++) {
        reducededge = sourceset[edgej] + " " + targetset[edgej];
        
        if (std::find(reducededgeset.begin(), reducededgeset.end(), reducededge) == reducededgeset.end()) {
            reducededgeset.push_back(reducededge);
            //cout << reducededge << endl;
            //edgefile << reducededge << "\n";
            for ( unsigned j = 0; j < eventind.size(); j++) {
                if (eventind[j] == sourceset[edgej]) nusource = eventnam[j];
                if (eventind[j] == targetset[edgej]) nutarget = eventnam[j];
            }
            //edgefile << nusource << "->" << nutarget << "\n";
            edgefile.push_back(nusource+"->"+nutarget);
        }
        
    }
    
    //edgefile.close();
    
    //return 0;
    return edgefile;
}
