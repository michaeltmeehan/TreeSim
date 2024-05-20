function isleaf(node::DataFrameRow)::Bool
    return node.left == 0 && node.right == 0
end

function isbinary(node::DataFrameRow)::Bool
    return node.right != 0
end

function isunary(node::DataFrameRow)::Bool
    return node.left != 0 && node.right == 0
end


function isroot(node::DataFrameRow)::Bool
    return node.host == 0
end

function binarize(tree::DataFrame)
    ctree = copy(tree)
    binary_tree = DataFrame()
    for row in eachrow(ctree)
        if isleaf(row)
            push!(binary_tree, row)
        elseif isbinary(row)
            left = row.left
            right = row.right
            while isunary(tree[left, :])
                left = tree[left, :left]
            end
            while isunary(tree[right, :])
                right = tree[right, :left]
            end
            row.left = left
            row.right = right
            push!(binary_tree, row)
        elseif isroot(row)
            left = row.left
            while isunary(tree[left, :])
                left = tree[left, :left]
            end
            row.left = left
            push!(binary_tree, row)
        end
    end
    new_ids = Dict(zip(binary_tree.id, 1:nrow(binary_tree)))
    new_ids[0] = 0
    binary_tree.id = [new_ids[id] for id in binary_tree.id]
    binary_tree.left = [new_ids[id] for id in binary_tree.left]
    binary_tree.right = [new_ids[id] for id in binary_tree.right]
    return binary_tree
end


function relabel!(wtree::DataFrame)::DataFrame
    new_ids = Dict{Int64, Int64}()
    for node in eachrow(wtree)
        if node.left == 0 && node.right == 0 # Leaf node
            if node.leaf_id == node.host # Sampled leaf
                new_ids[node.id] = CantorPair(node.id, node.host)
                node.id = new_ids[node.id]
            else # Transmission leaf
                new_ids[node.id] = CantorPair(0, node.leaf_id)
                node.id = new_ids[node.id]
            end
        else    # Internal node
            new_ids[node.id] = CantorPair(node.id, node.host)
            node.id = new_ids[node.id]
            node.left = new_ids[node.left]
            node.right = node.right == 0 ? 0 : new_ids[node.right]
        end
    end
    return wtree
end


"""
    normalize!(tree::DataFrame)::DataFrame
    Sequentially label nodes in `tree`.
"""
function normalize!(tree::DataFrame)::DataFrame
    new_labels = Dict(zip(tree.id, 1:nrow(tree)))
    new_labels[0] = 0
    tree.id = 1:nrow(tree)
    tree.left = [new_labels[node] for node in tree.left]
    tree.right = [new_labels[node] for node in tree.right]
    return tree
end


function join_wtrees(trees::Vector{DataFrame})::DataFrame
    tree = vcat(trees...)
    tree = tree[(tree.leaf_id .== 0 .|| tree.leaf_id .== tree.host), :] # Remove duplicate nodes
    sort!(tree, :t, rev=true)
    normalize!(tree)
    return tree
end