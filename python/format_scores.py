from matchms import Scores
from pandas import DataFrame
​
def create_long_table(data: DataFrame, value_id: str) -> DataFrame:
    """Convert the table from compact into long format.
    See DataFrame.melt(...).
    Args:
        data (DataFrame): The data table to convert.
        value_id (str): The name to assign to the added column through conversion to long format.
    Returns:
        DataFrame: Table in long format.
    """
    return data.transpose().melt(ignore_index=False, var_name='compound', value_name=value_id)
​
def join_df(x: DataFrame, y: DataFrame, on=[], how="inner") -> DataFrame:
    """Shortcut functions to join to dataframes on columns and index
    Args:
        x (DataFrame): Table X
        y (DataFrame): Table Y
        on (list, optional): Columns on which to join. Defaults to [].
        how (str, optional): Join method, see DataFrame.join(...). Defaults to "inner".
    Returns:
        DataFrame: Joined dataframe.
    """
    df_x = x.set_index([x.index] + on)
    df_y = y.set_index([y.index] + on)
    combined = df_x.join(df_y, how=how)
    return combined
​
def to_data_frame(scores) -> DataFrame:
    query_names = [spectra.get("id") for spectra in scores.queries]
    reference_names = [spectra.get("id") for spectra in scores.references]
​
    # Write scores to dataframe
    dataframe_scores = DataFrame(data=[entry["score"] for entry in scores.scores], index=reference_names, columns=query_names)
    dataframe_matches = DataFrame(data=[entry["matches"] for entry in scores.scores], index=reference_names, columns=query_names)
​
    scores_long = create_long_table(dataframe_scores, 'score')
    matches_long = create_long_table(dataframe_matches, 'matches')
​
    combined = join_df(matches_long, scores_long, on=['compound'], how='inner')
    return combined.reset_index().rename(columns={'level_0': 'query', 'compound': 'reference'})
