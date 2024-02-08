import pandas as pd
import argparse
import statsmodels.stats.multitest as smt


def choose_filters(column_names_str):
    filtering = True
    filters = []
    while filtering:
        another_filter = ""
        column = input("Enter the column name to filter: \nChoose from: " + column_names_str + "\n")
        column.lower()
        operator = input("Enter the operator (e.g. >, <): \n")
        value = input("Enter the value to filter: \n")
        filters.append((column, operator, value))
        while another_filter not in ["yes", "no", "y", "n"]:
            another_filter = input("Add another filter? (y/n)(yes/no)\n")
            if another_filter == "no" or another_filter == "n":
                filtering = False
    return filters


def filter(data_file, output_path):
    # Load the data from the provided file
    data = pd.read_csv(data_file, delimiter='\t')

    abbrev = {"chr": "Chromosome",
              "str": "Strand",
              "start": "Start",
              "end": "End",
              "id": "GeneID",
              "corcg": "Corrected_CG",
              "corchg": "Corrected_CHG",
              "corchh": "Corrected_CHH",
              "cor-noncg": "Corrected_non_CG",
              "len": "Gene length",
              "ex": "number of exons",
              "cy": "number of cytosines"}
    abbrev_column_names = [key + " (" + value + ")" for key, value in abbrev.items()]
    column_names_str = ',\n '.join(abbrev_column_names)

    done = False
    new = []
    not_found = []
    wrong_operator = []
    wrong_value = []
    filters = choose_filters(column_names_str)

    # Filter the data based on the provided filters
    filtered_data = data.copy()

    while not done:
        for column, operator, value in filters:
            try:
                filtered_data = filtered_data[
                    filtered_data[abbrev.get(column)].apply(lambda x: eval(f"x {operator} {value}"))]
            except KeyError:
                try:
                    filtered_data = filtered_data[filtered_data[column].apply(lambda x: eval(f"x {operator} {value}"))]
                except KeyError:
                    not_found.append(column)
                except TypeError and SyntaxError:
                    wrong_operator.append(operator)
                except NameError:
                    wrong_value.append(value)

            except (NameError, SyntaxError):
                try:
                    if column == "id" or column == "str":
                        filtered_data = filtered_data[filtered_data[abbrev.get(column)] == str(value)]
                except KeyError:
                    not_found.append(column)
                except TypeError and SyntaxError:
                    wrong_operator.append(operator)
                except NameError:
                    wrong_value.append(value)

            except TypeError:
                wrong_operator.append(operator)

        if not not_found and not wrong_operator and not wrong_value:
            done = True
        else:
            if not_found:
                print(f"Column(s) {','.join(not_found)} not found\n")
            elif wrong_operator:
                print(f"Operator(s) {','.join(wrong_operator)} not allowed for that column\n")
            elif wrong_value:
                print(f"Value(s) {','.join(wrong_value)} not allowed for that column\n")

            new_filter = input("Replace filter?(y/n)(yes/no)\n")

            if new_filter == "y" or new_filter == "yes":
                new += choose_filters(column_names_str)

            filters = new.copy()
            new.clear()
            not_found.clear()
            wrong_value.clear()
            wrong_operator.clear()

    filtered_data.to_csv(output_path, index=False, sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Command-line parser of Gene Methylation Classification")
    parser.add_argument("--data", help="Required Path to Data-File", required=True)
    parser.add_argument("--out", help="Required Path to Output-File", required=True)
    args = parser.parse_args()

    # Run filter method
    filter(args.data, args.out)
