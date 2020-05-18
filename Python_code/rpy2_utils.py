# external dependencies
import rpy2

python_str_to_R_str = rpy2.robjects.r('''
function(strg) {
    test <- strg
    test
}
''')