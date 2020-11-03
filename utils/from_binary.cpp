#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct VtxPair
{
    unsigned v;
    unsigned w;
};

auto read_word(std::ifstream & infile) -> unsigned
{
    unsigned char a, b;
    a = static_cast<unsigned char>(infile.get());
    b = static_cast<unsigned char>(infile.get());
    return unsigned(a) | (unsigned(b) << 8);
}

auto read_vf(const std::string & filename, bool unlabelled, bool undirected) -> void
{
    std::ifstream infile{ filename };
    if (! infile)
        throw "unable to open file";

    int size = read_word(infile);
    if (! infile)
        throw 0;

    std::vector<VtxPair> edges {};

    for (unsigned r = 0 ; r < size ; ++r) {
        read_word(infile);
    }

    if (! infile)
        throw 0;

	for (unsigned r = 0 ; r < size ; ++r) {
		int c_end = read_word(infile);
		if (! infile)
			throw 0;

		for (int c = 0 ; c < c_end ; ++c) {
			unsigned e = read_word(infile);

			if (e >= size)
				throw 1;

                        edges.push_back({r+1, e+1});
			read_word(infile);
		}
	}

    infile.peek();
    if (! infile.eof())
        throw 0;

    std::cout << "p edge " << size << " 0" << std::endl;

    for (auto e : edges)
        std::cout << "e " << e.v << " " << e.w << std::endl;
}

int main(int argc, char **argv) {
	read_vf(std::string {argv[1]}, true, true);
}
