struct similar {
    double fp;
    const char* patha;
    const char* pathb;
};

struct node {
    double fp;
    int patha;
    int pathb;

    struct node *next;
};

int append_to_list(double fp, int patha, int pathb);
int free_list();

int c_get_similarity_next(struct similar* a);

void c_setup( int N );
void c_add( int idx, const double a[], const char* b );

void c_double_to_float();
void c_teardown_floats();

void c_double_to_uint8();
void c_teardown_uint8();

void c_process();

void c_teardown();
