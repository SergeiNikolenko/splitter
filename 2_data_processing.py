import matplotlib.pyplot as plt
import pandas as pd


def load_data(file_path):
    return pd.read_csv(file_path, index_col=0)


def preprocess_data(info):
    # Заполнение пропущенных значений
    info["header"] = info["header"].fillna("")

    # Инициализация столбца для типов
    info["type"] = ""

    # Классификация по типам белков
    info.loc[info.header.str.contains("transport"), "type"] = "transport"
    info.loc[info.header.str.contains("hydrolase"), "type"] = "hydrolase"
    info.loc[info.header.str.contains("transferase"), "type"] = "transferase"
    info.loc[info.header.str.contains("transcription"), "type"] = "transcription"
    info.loc[info.header.str.contains("lyase"), "type"] = "lyase"
    info.loc[info.header.str.contains("oxidoreductase"), "type"] = "oxidoreductase"
    info.loc[info.header.str.contains("isomerase"), "type"] = "isomerase"
    info.loc[info.header.str.contains("ligase"), "type"] = "ligase"
    info.loc[info.header.str.contains("membrane"), "type"] = "membrane"
    info.loc[info.header.str.contains("viral"), "type"] = "viral"
    info.loc[info.header.str.contains("metal"), "type"] = "metal_containing"
    info.loc[info.header.str.contains("chaperone"), "type"] = "chaperone"

    # Остальные белки помечаются как "other"
    info.loc[info.type == "", "type"] = "other"

    return info


def save_data(info, output_path):
    info.to_csv(output_path)


def plot_type_distribution(info):
    # Статистика распределения типов белков
    type_counts = info.type.value_counts()
    type_counts.plot.pie(autopct="%1.1f%%", figsize=(7, 7))
    plt.title("Protein Type Distribution")
    plt.ylabel("")
    plt.show()


def main():
    # Загрузка данных
    input_file = "data/info.csv"
    output_file = "data/info_with_type.csv"

    info = load_data(input_file)

    # Предобработка данных
    info = preprocess_data(info)

    # Сохранение обработанных данных
    save_data(info, output_file)

    # Визуализация распределения типов белков
    plot_type_distribution(info)


if __name__ == "__main__":
    main()
