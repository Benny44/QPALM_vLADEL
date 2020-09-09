#include "ladel_types.h"
#include "ladel_constants.h"
#include "ladel_transpose.h"
#include "ladel_global.h"
#include "ladel_debug_print.h"

#define FIRST_LEAF 1
#define SUBSEQUENT_LEAF 2
#define NOT_A_LEAF -1

ladel_int ladel_least_common_ancestor(ladel_int subtree_root, ladel_int node, ladel_int *first_descendant,
                                        ladel_int *max_first_descendant, ladel_int *prev_leaf, 
                                        ladel_int *ancestor, ladel_int *node_type_of_leaf)
{
    if (subtree_root <= node || first_descendant[node] <= max_first_descendant[subtree_root])
    {
        *node_type_of_leaf = NOT_A_LEAF;
        return NONE;
    } else
    {
        max_first_descendant[subtree_root] = first_descendant[node];
        
        ladel_int last_leaf = prev_leaf[subtree_root];
        prev_leaf[subtree_root] = node;
        if (last_leaf == NONE)
        {
            *node_type_of_leaf = FIRST_LEAF;
            return node;
        } else
        {
            *node_type_of_leaf = SUBSEQUENT_LEAF;
            
            ladel_int lca; /*least common ancestor*/
            for (lca = last_leaf; lca != ancestor[lca]; lca = ancestor[lca]);
            
            ladel_int last_path_node, ancestor_of_last_path_node;
            for (last_path_node = last_leaf; last_path_node != lca; last_path_node = ancestor_of_last_path_node)
            {
                ancestor_of_last_path_node = ancestor[last_path_node];
                ancestor[last_path_node] = lca;
            }
            return lca;
        }  
    }
    
}

ladel_int ladel_col_counts(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_work *work)
{
    if (!M || !sym || !sym->etree || !sym->postorder || !sym->col_counts || !work) return FAIL;
    
    ladel_int *etree = sym->etree, *postorder = sym->postorder, *col_counts = sym->col_counts;
    ladel_int *prev_leaf = work->array_int_ncol1;
    ladel_int *first_descendant = work->array_int_ncol2;
    ladel_int *max_first_descendant = work->array_int_ncol3;
    ladel_int *ancestor = work->array_int_ncol4;
    ladel_int ncol = M->ncol, index, node, post_node, subtree_root, lca, type_of_leaf = NOT_A_LEAF;
    
    ladel_sparse_matrix *M_lower;
    if (M->symmetry == UPPER) M_lower = ladel_transpose(M, FALSE, work);
    else if (M->symmetry == LOWER) M_lower = M;
    else return FAIL;

    if (!M_lower) return FAIL;

    for (index = 0; index < ncol; index++) prev_leaf[index] = NONE;
    for (index = 0; index < ncol; index++) first_descendant[index] = NONE;
    for (index = 0; index < ncol; index++) max_first_descendant[index] = NONE;
    for (index = 0; index < ncol; index++) ancestor[index] = index;

    for (node = 0; node < ncol; node++)
    {
        post_node = postorder[node];
        if (first_descendant[post_node] == NONE) 
            col_counts[post_node] = 1;
        else 
            col_counts[post_node] = 0;
        for (; post_node != NONE && first_descendant[post_node] == NONE; post_node = etree[post_node]) 
            first_descendant[post_node] = node;
    }
    for (node = 0; node < ncol; node++)
    {
        post_node = postorder[node];
        if (!IS_ROOT(post_node, etree)) col_counts[etree[post_node]]--;
        LADEL_FOR(index, M_lower, post_node)
        {
            subtree_root = M_lower->i[index];
            lca = ladel_least_common_ancestor(subtree_root, post_node, first_descendant, max_first_descendant, prev_leaf, ancestor, &type_of_leaf);
            if (type_of_leaf != NOT_A_LEAF) col_counts[post_node]++;
            if (type_of_leaf == SUBSEQUENT_LEAF) col_counts[lca]--; /*correct for duplicates*/
        }
        if (!IS_ROOT(post_node, etree)) ancestor[post_node] = etree[post_node];
    }
    for (node = 0; node < ncol; node++)
        if (!IS_ROOT(node, etree)) col_counts[etree[node]] += col_counts[node];

    for (node = 1; node < ncol; node++)
    {
        col_counts[node] += --col_counts[node-1];
    }
    col_counts[ncol-1]--;
        
    if (M->symmetry == UPPER) ladel_sparse_free(M_lower);
    return col_counts[ncol-1];
}
