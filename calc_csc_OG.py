import pandas as pd

# Process Orthogroups.tsv
def read_orthogroups_results(filepath):
    dtype = {0: str}  # Assuming the first column is the only one that should be string
    
    try:
        df = pd.read_csv(filepath, sep='\t', dtype=dtype, low_memory=False)
        
        # Exclude rows where 'Orthogroup' column contains the word 'Orthogroup'
        df = df[~df['Orthogroup'].str.contains("Orthogroup", na=False)]
        
        return df
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        return None
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{filepath}' is empty.")
        return None
    except pd.errors.ParserError:
        print(f"Error: There was an error parsing the file '{filepath}'.")
        return None

# Calculate percentage of each category in totals, ensuring sum to 100%
def calculate_percentage_totals(df_totals):
    relevant_totals = df_totals.loc[df_totals["Category"].isin(["core", "shell", "cloud"])]
    column_sums = relevant_totals.select_dtypes(include=[int, float]).sum(axis=0)
    percentage_totals = relevant_totals.select_dtypes(include=[int, float]).div(column_sums, axis=1) * 100
    percentage_totals.insert(0, "Category", relevant_totals["Category"].values)
    return percentage_totals


def calculate_category_percentages(gene_counts):
    # Exclude specified columns
    exclude_columns = ['Orthogroup', 'Zero_Count', 'category', 'total']
    relevant_columns = [col for col in gene_counts.columns if col not in exclude_columns]

    # Initialize a DataFrame to hold percentage results
    percentage_results = pd.DataFrame()

    for col in relevant_columns:
        # Group by category and sum the values for each category
        category_sums = gene_counts.groupby('category')[col].sum()
        
        # Calculate the total for the column
        total_sum = category_sums.sum()
        
        # Calculate percentages
        percentages = (category_sums / total_sum) * 100
        
        # Add results to the percentage_results DataFrame
        percentage_results[col] = percentages

    # Reset index to have categories as a column
    percentage_results = percentage_results.reset_index()
    percentage_results.columns = ['Category'] + relevant_columns  # Rename columns

    return percentage_results


def main():
    orthogroups_filepath = "Orthogroups.complete.tsv"

    # Process Orthogroups.tsv (OG)
    df_OG = read_orthogroups_results(orthogroups_filepath)

    if df_OG is not None:
        # Display the first few rows of the DataFrame
        print("Orthogroups Data:")
        print(df_OG.head())
        
        # Select the Orthogroup column and all other columns dynamically
        gene_counts_columns = df_OG.columns[1:]  # Get all columns except the first one
        gene_counts = df_OG[['Orthogroup']].copy()  # Start with the Orthogroup column
        
        # Calculate gene counts for all other columns
        for col in gene_counts_columns:
            gene_counts[col] = df_OG[col].map(lambda cell: 0 if pd.isna(cell) or cell == "" else len(cell.split(",")))

        # Count the number of zeros per row (excluding header)
        gene_counts['Zero_Count'] = gene_counts[gene_counts_columns].apply(lambda row: (row == 0).sum(), axis=1)

        # Determine the number of relevant columns for the category condition
        num_relevant_columns = len(gene_counts_columns) - 1  # Excluding 'Orthogroup'

        # Add the category column based on the conditions specified
        gene_counts['category'] = gene_counts.apply(
            lambda row: 'core' if row['Zero_Count'] == 0 
            else 'cloud' if row['Zero_Count'] == num_relevant_columns 
            else 'shell',  # Assign 'shell' if neither condition is met
            axis=1
        )

        # Add a total column by summing across the gene counts
        gene_counts['total'] = gene_counts[gene_counts_columns].sum(axis=1)

        # Create a new DataFrame for each category summary
        core_total = gene_counts[gene_counts['category'] == 'core']['total'].sum()
        shell_total = gene_counts[gene_counts['category'] == 'shell']['total'].sum()
        cloud_total = gene_counts[gene_counts['category'] == 'cloud']['total'].sum()

        # Prepare summary DataFrame with categories
        summary_data = {
            'Category': ['core', 'shell', 'cloud'],
            'Total': [core_total, shell_total, cloud_total]
        }
        summary_df = pd.DataFrame(summary_data)

        # Calculate percentage of each category
        percentage_summary = calculate_percentage_totals(summary_df)

        # New DataFrame to sum totals by Zero_Count
        max_zero_count = gene_counts['Zero_Count'].max()
        zero_count_summary = pd.DataFrame({'Zero_Count': range(max_zero_count + 1), 'Total': 0})

        for zero_count in range(max_zero_count + 1):
            total_sum = gene_counts[gene_counts['Zero_Count'] == zero_count]['total'].sum()
            zero_count_summary.loc[zero_count_summary['Zero_Count'] == zero_count, 'Total'] = total_sum

        # Add category column to zero_count_summary
        zero_count_summary['category'] = zero_count_summary['Zero_Count'].apply(
            lambda x: 'core' if x == 0 else 'cloud' if x == max_zero_count else 'shell'
        )

        # Calculate percentage of core, shell, and cloud per column
        category_percentage_df = calculate_category_percentages(gene_counts)

        # Display results
        print("\nGene Counts with Zero Count, Category, and Total:")
        print(gene_counts)

        print("\nSummary DataFrame with Totals:")
        print(summary_df)

        print("\nPercentage Summary DataFrame:")
        print(percentage_summary)

        print("\nZero Count Summary DataFrame:")
        print(zero_count_summary)

        print("\nCategory Percentage DataFrame:")
        print(category_percentage_df)

        # Save the DataFrames to new CSV files
        gene_counts.to_csv("Gene_Counts.csv", index=False)
        summary_df.to_csv("Summary.csv", index=False)
        percentage_summary.to_csv("Percentage_Summary.csv", index=False)
        zero_count_summary.to_csv("Zero_Count_Summary.csv", index=False)
        category_percentage_df.to_csv("Category_Percentage_Summary.csv", index=False)

if __name__ == "__main__":
    main()
