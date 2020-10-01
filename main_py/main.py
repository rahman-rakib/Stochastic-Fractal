
def test_number_length():
    from main_py.stochastic_fragmentation import NumberLength
    sfl = NumberLength(alpha=2, probability=0.75)
    sfl.log(True)
    sfl.run(100, 20, 10)
    sfl.view()
    print(sfl.get_signature())
    print(sfl.get_filename())
    sfl.run_ensemble(10, 1000, 500, 10)
    pass


def test_true_length():
    from main_py.stochastic_fragmentation import TrueLengths
    sfl = TrueLengths(alpha=2, probability=0.75)
    sfl.log(True)
    sfl.run(100)
    sfl.view()
    print(sfl.get_signature())
    print(sfl.get_filename())
    sfl.run_ensemble(10, 1000)
    pass


def test_moment():
    from main_py.stochastic_fragmentation import Moment
    sfl = Moment(alpha=2, probability=0.75, exponent=0.7)
    sfl.log(True)
    a = sfl.run(1000, 500, 10)
    print("moment a ", a)
    sfl.view()
    print(sfl.get_signature())
    print(sfl.get_filename())
    sfl.run_ensemble(10, 1000, 500, 10)
    pass


def test_df_determination():
    from main_py import df_determination
    df_determination.test()


def test_analytic():
    from  main_py import analytic_solution
    # analytic_solution.test()
    analytic_solution.test_2()


def main():
    # test_number_length()
    # test_true_length()
    # test_moment()
    # test_df_determination()
    test_analytic()
    pass


if __name__ == '__main__':
    main()