import process

for lmax in range(1, 41):
    file_path = f'coefficients/degree_1_to_{lmax}.coef'
    process.convert_to_sh(file_path, lmax)