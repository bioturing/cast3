#include "generate.h"
#include <assert.h>

int8_t ascii_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

double prob_barcode_err[16][4] = {
    { 0.00243200183210607,   0.00265226825720049,   0.00238252487266703,   0.00247859241604291  },
    { 9.84518532280806e-5,   0.000105418767099898,  0.00012024540587624,   0.000149312364560738 },
    { 6.43035178977319e-5,   8.69887571314547e-05,  7.63701208403244e-5,   7.50233544318644e-5  },
    { 7.23176133316617e-5,   8.47922299523897e-05,  7.09183994956776e-5,   9.2661673292949e-5   },
    { 6.0728820107815e-5,    8.68149497299812e-05,  6.71073268524629e-5,   0.000104612297544355 },
    { 5.75798861973855e-5,   6.53444346828966e-05,  6.45765176008781e-5,   9.95650559011039e-5  },
    { 0.000113841438435535,  0.000140194463619207,  0.000193938780289159,  0.000151948068667371 },
    { 0.000337649068858504,  0.000249891993056765,  0.000186216398571565,  0.000232127896750806 },
    { 0.000208871197414734,  0.000252154220076128,  0.000218352053830969,  0.000268515970157629 },
    { 0.000220690984209116,  0.000171472165707382,  0.000161881515059254,  0.000244160685961411 },
    { 0.000217392257875513,  0.000222952327734294,  0.000193974762597782,  0.000228887240080449 },
    { 0.000200543541862804,  0.000196424451019658,  0.000162242944498548,  0.000226064717112528 },
    { 0.000204301845807882,  0.00021033089044356,   0.000204731464934722,  0.000255853410228468 },
    { 0.000201081348868472,  0.000242918171866774,  0.000236683915627231,  0.000240610244102754 },
    { 0.000213428902961876,  0.00025577076336335,   0.000205504040440625,  0.000223785784021304 },
    { 0.000376832758819433,  0.000361701876772859,  0.000321297197848344,  0.000475347661405257 }
};

double prob_barcode_mut[16][4][5] = {
    /* 0 */
    { {0, 0.0150060806230109, 0.0372176929088549, 0.0456575461284543, 1},
      {0.018543294755385, 0, 0.038966807605046, 0.0544988281077176, 1},
      {0.0199607426916924, 0.0283670224586378, 0, 0.0606673270456195, 1},
      {0.0114106115349607, 0.0240909626615095, 0.066435695729253, 0, 1} },

    /* 1 */
    { {0, 0.248391841609553, 0.452225152780146, 0.563918024508456, 1},
      {0.271759809161439, 0, 0.383139952595116, 0.596635639908786, 1},
      {0.253201640746953, 0.356231919488245, 0, 0.612659080115662, 1},
      {0.0867939772665535, 0.315210731658511, 0.681413839419953, 0, 1} },

    /* 2 */
    { {0, 0.35875670578243, 0.646852731962304, 0.73066954778514, 1},
      {0.343760589219978, 0, 0.479971081907819, 0.793895808076938, 1},
      {0.285338081360802, 0.396841355796019, 0, 0.763238404334643, 1},
      {0.105581623902801, 0.316107883346555, 0.785941934496185, 0, 1} },

    /* 3 */
    { {0, 0.381603629518156, 0.720009551363568, 0.809796811399441, 1},
      {0.310711927790844, 0, 0.494214799315343, 0.833901512819323, 1},
      {0.304626183783441, 0.425403239506486, 0, 0.811649040628603, 1},
      {0.115669720635526, 0.338093421976441, 0.854613457687941, 0, 1} },

    /* 4 */
    { {0, 0.422569262592116, 0.8786178694049, 0.993668893861723, 1},
      {0.418755753334043, 0, 0.699391707797705, 0.994838536583109, 1},
      {0.416612808249673, 0.598519013218052, 0, 0.994178503549868, 1},
      {0.171887212260899, 0.408418695997334, 0.995999176955996, 0, 1} },

    /* 5 */
    { {0, 0.455596441354606, 0.890908380782869, 0.999912121879978, 1},
      {0.413720923086378, 0, 0.692969916725563, 0.999904126847556, 1},
      {0.385448512090586, 0.537184675888728, 0, 0.999902986763912, 1},
      {0.143105140923864, 0.376153257273678, 0.99994272535355, 0, 1} },

    /* 6 */
    { {0, 0.406327122386922, 0.826869104480771, 0.997881317068697, 1},
      {0.374638212546548, 0, 0.705784016041249, 0.998169578917216, 1},
      {0.310655760448629, 0.673039164730911, 0, 0.998758410649673, 1},
      {0.238939465779065, 0.481877157622821, 0.998489829208511, 0, 1} },

    /* 7 */
    { {0, 0.195755343642281, 0.597418219800039, 0.9999027098749, 1},
      {0.306706272884354, 0, 0.623375596535977, 0.999887185275846, 1},
      {0.425573183339738, 0.644636465230306, 0, 0.999868880685305, 1},
      {0.268149141057318, 0.583011172880877, 0.999894468159314, 0, 1} },

    /* 8 */
    { {0, 0.518582892018808, 0.784930884223576, 0.999240164825784, 1},
      {0.27197526200349, 0, 0.677298530767374, 0.999459461525574, 1},
      {0.298843598520712, 0.566920693914455, 0, 0.999250351099056, 1},
      {0.206940109327211, 0.453842829047238, 0.999488510173712, 0, 1} },

    /* 9 */
    { {0, 0.310375812444249, 0.711274446114916, 0.995279006975584, 1},
      {0.33198496434301, 0, 0.682154642434746, 0.994291602749511, 1},
      {0.326298575650106, 0.623180862865886, 0, 0.993958370958107, 1},
      {0.245572771549946, 0.509869938136896, 0.995911758726658, 0, 1} },

    /* 10 */
    { {0, 0.328380850761973, 0.787038692763321, 0.999824876157147, 1},
      {0.314533542852326, 0, 0.761671785760346, 0.999832485674823, 1},
      {0.394375376537983, 0.67886273624847, 0, 0.999770609437724, 1},
      {0.215190322042968, 0.459559438903071, 0.999832267459712, 0, 1} },

    /* 11 */
    { {0, 0.351088558521034, 0.789439220146743, 0.998694370574477, 1},
      {0.313613251580904, 0, 0.711908553951083, 0.998681302505116, 1},
      {0.359368440646015, 0.617058478744723, 0, 0.998362884254382, 1},
      {0.209646159907881, 0.468978015961584, 0.998889020344748, 0, 1} },

    /* 12 */
    { {0, 0.393658697481948, 0.792198525518021, 0.999928449922612, 1},
      {0.410475580222229, 0, 0.722234908570464, 0.999931264601499, 1},
      {0.374465971655011, 0.658098282400073, 0, 0.99992428468644, 1},
      {0.253478004683067, 0.534648200183392, 0.99995165622933, 0, 1} },

    /* 13 */
    { {0, 0.356711759499212, 0.728725231428964, 0.999992410855512, 1},
      {0.382426264718291, 0, 0.663486603263781, 0.999997354906714, 1},
      {0.35486390882733, 0.602597758689963, 0, 0.999991516359174, 1},
      {0.232284719632891, 0.495108208909138, 0.999995994302563, 0, 1} },

    /* 14 */
    { {0, 0.319898694544483, 0.731679405112689, 1, 1},
      {0.37096552352182, 0, 0.706085746854284, 1, 1},
      {0.390778933035154, 0.642144404001968, 0, 1, 1},
      {0.273646850490209, 0.52192922481363, 1, 0, 1} },

    /* 15 */
    { {0, 0.2932143024609, 0.710902898260362, 0.999975702185524, 1},
      {0.354878944582217, 0, 0.677907230816238, 0.999974463693145, 1},
      {0.383118046981612, 0.670801741155919, 0, 0.999968252627095, 1},
      {0.258755909167866, 0.495252896670793, 0.999974824047124, 0, 1} }
};

double prob_seq_mut[4][5] = {
    {0, 0.381603629518156, 0.720009551363568, 0.809796811399441, 1},
    {0.310711927790844, 0, 0.494214799315343, 0.833901512819323, 1},
    {0.304626183783441, 0.425403239506486, 0, 0.811649040628603, 1},
    {0.115669720635526, 0.338093421976441, 0.854613457687941, 0, 1}
};

std::random_device rd;
std::mt19937 generator(rd());
std::string ascii_char = "ACGTN";

generate_t::generate_t(const char *bx_path, const char *hap1_path, const char *bed1_path,
                       const char *hap2_path, const char *bed2_path, const char *fai_path)
: hap1(hap1_path, bed1_path, fai_path, 1), hap2(hap2_path, bed2_path, fai_path, 2)
{
    fprintf(stderr, "Loading barcode ... ");
    std::ifstream fi(bx_path);
    if (!fi.is_open()) {
        perror("\nCould not open barcode file");
        exit(EXIT_FAILURE);
    }
    std::string bx;
    while (fi >> bx)
        barcode_lst.push_back(bx);
    fi.close();
    fprintf(stderr, "done\n");

    fi_mlc.open("molecule.inf");
}

molecule_t generate_t::generate_molecule(int read_len)
{
    std::gamma_distribution<double> gamma_dis(GAMMA_DIS_A, GAMMA_DIS_B);
    std::poisson_distribution<int> poisson_dis(MEAN_MLC_COVERAGE);
    std::uniform_int_distribution<int> uniform_01(0, 1);
    int mole_len, coverage, n_read_pair;

    mole_len = (int)gamma_dis(generator);
    mole_len = std::min(MAX_MLC_LEN, mole_len);
    mole_len = std::max(MIN_MLC_LEN, mole_len);
    coverage = poisson_dis(generator);
    n_read_pair = std::max(mole_len * coverage / 100 / (read_len * 2 - BARCODE_SZ), 1);
    molecule_t ret(n_read_pair, mole_len);

    /* get pos */
    ret.hap = uniform_01(generator);
    std::uniform_int_distribution<int64_t> uniform_dis(0, ret.hap ?
                                                       hap1.get_chr_pos_back():
                                                       hap2.get_chr_pos_back());
    std::pair<int, int> rp;
    while (1) {
        rp = ret.hap ? hap1.check_pos(uniform_dis(generator), mole_len) :
                       hap2.check_pos(uniform_dis(generator), mole_len);
        if (rp.first != -1)
            break;
    }
    ret.tid = rp.first;
    ret.start_pos = rp.second;
    /* ret.start_pos is base-1, chr is base-0 */
    assert(ret.start_pos - 1 + mole_len <= (ret.hap ? hap1.get_chr_size(ret.tid) : 
                                                      hap2.get_chr_size(ret.tid)));
    return ret;
}

std::vector<molecule_t> generate_t::generate_barcode(int read_len, int &total_read,
                                                     int mean_mlc_per_bx)
{
    std::poisson_distribution<int> poisson_dis(mean_mlc_per_bx);
    int n_mole = poisson_dis(generator);
    std::vector<molecule_t> ret;
    while (n_mole--) {
        ret.push_back(generate_molecule(read_len));
        total_read -= ret.back().n_read_pair;
    }
    return ret;
}

std::string gen_barcode_error(const std::string &barcode)
{
    std::string ret = barcode;
    std::uniform_real_distribution<double> uniform_dis(0.0, 1.0);
    int i, j;

    for (i = 0; i < (int)ret.size(); ++i) {
        int8_t c = ascii_table[(int)ret[i]];
        assert(c < 4);
        if (uniform_dis(generator) > prob_barcode_err[i][c])
            continue;
        for (j = 0; j < 5; ++j) {
            if (uniform_dis(generator) > prob_barcode_mut[i][c][j])
                continue;
            ret[i] = ascii_char[j];
            break;
        }
    }
    return ret;
}

void gen_read_error(std::string &seq1, std::string &seq2)
{
    std::uniform_real_distribution<double> uniform_dis(0.0, 1.0);
    std::uniform_int_distribution<int> uniform_03(0, 3);
    int i, j;

    for (i = 0; i < (int)seq1.size(); ++i) {
        seq1[i] = toupper(seq1[i]);
        if (ascii_table[(int)seq1[i]] == 4)
            seq1[i] = ascii_char[(uniform_03(generator))];
        if (uniform_dis(generator) > READ1_ERR)
            continue;
        int8_t c = ascii_table[(int)seq1[i]];
        assert(c < 4);
        for (j = 0; j < 5; ++j) {
            if (uniform_dis(generator) > prob_seq_mut[c][j])
                continue;
            seq1[i] = ascii_char[j];
            break;
        }
    }

    for (i = 0; i < (int)seq2.size(); ++i) {
        seq2[i] = toupper(seq2[i]);
        if (ascii_table[(int)seq2[i]] == 4)
            seq2[i] = ascii_char[(uniform_03(generator))];
        if (uniform_dis(generator) > READ2_ERR)
            continue;
        int8_t c = ascii_table[(int)seq2[i]];
        assert(c < 4);
        for (j = 0; j < 5; ++j) {
            if (uniform_dis(generator) > prob_seq_mut[c][j])
                continue;
            seq2[i] = ascii_char[j];
            break;
        }
    }
}

void generate_t::generate_read(int read_len, int total_read, int mean_mlc_per_bx)
{
    std::normal_distribution<double> normal_dis(NORMAL_MEAN, NORMAL_STDDEV);
    int i, j, k, mlc_id = 0;
    fi_mlc << "#mlc_id\thaplotype\tref\tstart_pos\tlen\tn_read_pair\n";

    for (i = 0; ; ++i) {
        if (i == (int)barcode_lst.size()) {
            fprintf(stderr, "Total barcodes consume is exceed limit!\n"
                    "Please decrease # million reads pairs in total to simulated!\n");
            exit(EXIT_FAILURE);
        }

        auto ret = generate_barcode(read_len, total_read, mean_mlc_per_bx);
        std::string seq1, seq2;
        int pos, isize;

        for (auto &mole : ret) {
            genome_t &hap = mole.hap ? hap1 : hap2;

            /* output molecule data */
            fi_mlc << mlc_id << "\t" << (mole.hap ? "1" : "2") << "\t"
                   << hap.get_ref_name(mole.tid) << "\t" << mole.start_pos << "\t"
                   << mole.len << "\t" << mole.n_read_pair << "\n";

            for (j = 0; j < mole.n_read_pair; ++j) {
                std::uniform_int_distribution<int> uniform_dis(0, mole.len - 1);
                std::string barcode = gen_barcode_error(barcode_lst[i]);
                while (1) {
                    pos = uniform_dis(generator);
                    isize = (int)normal_dis(generator);
                    if (pos + isize + read_len <= mole.len)
                        break;
                }

                /* ret.start_pos is base-1, chr is base-0 */
                seq1 = hap.get_sub_chr(mole.tid, mole.start_pos - 1 + pos, read_len - BARCODE_SZ);
                seq2 = hap.get_sub_chr(mole.tid, mole.start_pos - 1 + pos + isize, read_len);
                gen_read_error(seq1, seq2);

                /* output read */
                std::string name = hap.convert(mole.tid, mole.start_pos + pos, read_len - BARCODE_SZ,
                                               mole.start_pos + pos + isize, read_len);
                name = std::to_string(mlc_id) + "|" + name;

                printf("@%s/1\tBX:%s\tQB:", name.c_str(), barcode.c_str());
                for (k = 0; k < BARCODE_SZ; ++k)
                    printf("#");
                printf("\n");
                printf("%s\n+\n", seq1.c_str());
                for (k = 0; k < read_len; ++k)
                    printf("#");
                printf("\n");

                printf("@%s/2\tBX:%s\tQB:", name.c_str(), barcode.c_str());
                for (k = 0; k < BARCODE_SZ; ++k)
                    printf("#");
                printf("\n");
                printf("%s\n+\n", seq2.c_str());
                for (k = 0; k < read_len; ++k)
                    printf("#");
                printf("\n");
            }
            ++mlc_id;
        }

        if (total_read <= 0)
            break;
    }
}
