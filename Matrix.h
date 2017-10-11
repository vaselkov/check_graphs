
class Matrix{
    private:
        int* data;
        unsigned msize;
    public:
        Matrix();
        Matrix(unsigned);
        Matrix(const Matrix&);
        ~Matrix();
        int* operator[](unsigned) const;
        //ifstream& operator >> (ifstream& in);
        unsigned size() const;
        void show(std::ofstream&);
        void show();
        void clear();
        void create(unsigned);
        void create(unsigned, int);
};

