
def test_number_length():
    from main_py.stochastic_fragmentation import NumberLength
    sfl = NumberLength(alpha=2, probability=0.75)
    sfl.log(True)
    sfl.run(100, 20, 10)
    sfl.view()
    pass

def main():
    test_number_length()
    pass


if __name__ == '__main__':
    main()