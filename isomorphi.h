#include <vector>

bool isomorphi(const Matrix &X, const Matrix &Y,
        bool (*Compare)(const Matrix&, const Matrix&, const std::vector<int>&),
        void (*getDeg)(const Matrix&, std::vector<int>&));
bool isomorph(const Matrix &X, const Matrix &Y,
        bool (*Compare)(const Matrix&, const Matrix&, const std::vector<int>&));
