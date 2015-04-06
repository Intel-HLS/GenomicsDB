#ifndef PARSE_COMMAND_LINE_H
#define PARSE_COMMAND_LINE_H

#include<fstream>
class CommandLineOpts
{
  public:
    CommandLineOpts()
    {
      m_do_scan = false;
      m_is_input_csv_sorted = false;
      m_workspace = 0;
      m_csv_filename = 0;
      m_array_name = 0;
      m_num_samples = 0ull;
      m_position = 0ull;
      m_temp_space = "";
    }
    bool m_do_scan;
    bool m_is_input_csv_sorted;
    char* m_workspace;
    char* m_csv_filename;
    char* m_array_name;
    std::ofstream m_output_fstream;
    std::ifstream m_positions_list;
    uint64_t m_num_samples;
    uint64_t m_position;
    std::string m_temp_space;
};

void parse_command_line(int argc, char** argv, CommandLineOpts& cl);

#endif
