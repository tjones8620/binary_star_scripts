from texttable import Texttable
import numpy as np
import latextable
import argparse

parser = argparse.ArgumentParser(description='Convert texttable to latex table.')
parser.add_argument('path', type=str, help='Full Path to csv file.')
parser.add_argument('--caption', type=str, help='Caption for latex table.', default='Caption')
parser.add_argument('--label', type=str, help='Label for latex table.', default='Label')
args = parser.parse_args()

def main():
        data = np.loadtxt(args.path, dtype=str, delimiter=',').tolist()

        table = Texttable()
        table.set_cols_align(["c"] * len(data[0]))
        table.set_deco(Texttable.HEADER | Texttable.BORDER)
        table.add_rows(data)

        print('\nTexttable Table:')
        print(table.draw())

        print('\nTexttable Latex:')
        print(latextable.draw_latex(table, caption=args.caption, label=args.label))

if __name__ == '__main__':
        main()

