import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import math

st.set_page_config(layout="wide")
st.title("Kinetics Lab 1 Calculator")

a = 1.38e5  # L / g*cm (absorptivity const. of red40)
b = 1  # cm

graph_order = 0

bleach_conc0 = 1
bleach_order = 1


def process_data(df):
    # steps:
    # 1. get concentration
    # 2. graph for order 0..2
    # 3. determine order and calculate k
    st.markdown(
        r"To find concentration, we can use the formula $A = abc$. We know $A$, $a$ and $b$, so we can find $c$ with $c = \frac{A}{ab}$. Doing this for each data point yields the following data:"
    )
    abs_col = df["absorptivity"]
    concs = [x / (a * b) for x in abs_col]

    # what order is the rxn with respect to red40
    most_linear = -1
    most_linear_r2 = 0
    slope = 0

    cols = st.columns(3)
    fns = [
        lambda x: x,
        lambda x: [math.log(n) for n in x],
        lambda x: [1 / n for n in x],
    ]
    for idx, (col_ctx, fn) in enumerate(zip(cols, fns)):
        if graph_order != -1 and idx != graph_order:
            continue
        with col_ctx:
            fig, ax = plt.subplots()
            ax.plot(df["time"], fn(concs))
            ax.set_xlabel("Time")
            ax.set_ylabel("Concentration of red40 (g/L)")
            ax.set_title(f"Order {idx}")
            st.pyplot(fig)

            # linear regress to determine if the graph is straight
            m, _, r, _, _ = stats.linregress(df["time"], fn(concs))
            if r * r > most_linear_r2:
                most_linear = idx
                most_linear_r2 = r * r
                slope = m

    # check if the r squared is "linear"
    st.markdown(f"r^2^: {most_linear_r2}")
    if most_linear_r2 < 0.9:
        st.markdown(
            "Note: r^2^ is less than 0.9 for the most linear plot. Data may be noisy."
        )

    bleach_const = bleach_conc0**bleach_order

    st.markdown(f"Order {most_linear} graph was found to be most linear.")

    st.markdown(f"The rate law for this reaction is:")
    st.latex(f"r=k[bleach]^m[red40]^{most_linear}")
    st.markdown(
        "Assuming we have a large excess of bleach, we can determine that the concentration of bleach remains relatively constant to red40."
    )
    st.markdown(
        f"This means that we can simplify the rate law to $r=k'\\cdot[red40]^{most_linear}$, where k' is the pseudo-rate constant whose value is equal to $k\\cdot[bleach]^m_0 = k\\cdot{bleach_const}$"
    )

    if most_linear == 2:
        kprime = slope
    else:
        kprime = -slope

    k = kprime / bleach_const
    st.markdown(
        f"From the graph, we can determine that k' = {kprime:.3f}. We can derive that $k = \\frac{{k'}}{{[bleach]^m_0}} = \\frac{{{kprime:.3f}}}{{{bleach_const}}} = {k}$. Our final answers:"
    )
    st.latex(f"k = {k:.4f}")
    st.latex(f"r = {k:.4f}\\cdot[bleach]^{bleach_order}[red40]^{most_linear}")


def main():
    global bleach_order, bleach_conc0, graph_order
    bleach_order = st.number_input(
        "Order of reaction in respect to bleach (if given): ", value=1
    )
    graph_order = st.number_input("Graph order (-1 for all):", min_value=-1, value=-1, max_value=2)

    data = st.file_uploader("Add absorptivity data here:", type=["csv"])
    if data is None:
        return

    # convert to lowercase to avoid captilization issues
    df = pd.read_csv(data)
    df.columns = df.columns.str.lower()
    st.dataframe(df)

    # check file integrity
    required_cols = ["time", "absorptivity"]
    missing = [col for col in required_cols if col not in df.columns]
    if len(missing) > 0:
        st.error(f"Missing columns: {missing}")
        st.stop()

    process_data(df)


main()
