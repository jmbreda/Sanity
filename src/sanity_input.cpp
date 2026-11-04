#include "sanity_input.hpp"

#include <map>

#include "logging.hpp"

CountIndex load_count_index(const SanityOptions &options)
{
    CountIndex counts;
    logging_debug("Reading input file to get gene and cell counts");
    if (options.in_file_extension == "mtx")
    {
        std::map<int, int> row_to_gene; // Maps each MTX row to the index of its gene among the kept genes (-1 if the row has no counts).
        index_mtx_counts(options.in_file, counts.num_input_rows, counts.num_genes, counts.num_cells, row_to_gene,
                         counts.mtx_rows_g, counts.cell_totals_c, counts.gene_totals_g, options.num_threads);
        counts.gene_names_g = Read_GeneNames(options.gene_name_file, row_to_gene, counts.num_genes);
        if (counts.gene_names_g.empty())
        {
            logging_debug("Warning: gene names file not found or empty, generating default gene names.");
            for (int g = 0; g < counts.num_genes; ++g)
                counts.gene_names_g.push_back("Gene_" + std::to_string(g + 1));
        }
        counts.cell_names_c = Read_CellNames(options.cell_name_file);
        if (counts.cell_names_c.empty())
        {
            logging_debug("Warning: cell names file not found or empty, generating default cell names.");
            for (int c = 0; c < counts.num_cells; ++c)
                counts.cell_names_c.push_back("Cell_" + std::to_string(c + 1));
        }
    }
    else
    {
        index_tsv_counts(options.in_file, counts.num_input_rows, counts.num_genes, counts.num_cells, counts.tsv_offsets_g,
                         counts.cell_totals_c, counts.gene_totals_g, counts.cell_names_c, counts.gene_names_g, options.num_threads);
    }
    logging_debug("There were " + std::to_string(counts.num_input_rows) + " rows");

    std::string genes_message = "The input contained " + std::to_string(counts.num_input_rows) + " rows for " + std::to_string(counts.num_cells) +
                                " cells, and " + std::to_string(counts.num_genes) + " of these genes had at least one count in the dataset.";
    if (counts.num_input_rows != counts.num_genes)
        genes_message += " " + std::to_string(counts.num_input_rows - counts.num_genes) + " genes were thus removed because no counts were detected.";
    logging_debug(genes_message);
    return counts;
}
