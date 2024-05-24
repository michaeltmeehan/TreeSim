@recipe function f(phylo::Phylogeny)

    timespan = extrema([node.data["t"] for node in traversal(phylo.tree, inorder)])
    xlims := timespan
    xscale := :identity
    legend := nothing
    xlabel := "Time"
    yaxis := false

    @series begin
        linewidth := 5
        linecolor := [node.data["host"] for node in traversal(phylo.tree, preorder)]
        # markersize := 10
        # markercolor := [node.id for node in traversal(phylo.tree, preorder)]
        # markerstrokewidth := 0
        seriestype := :path
        phylo.tree
    end

    # @series begin
    #     markersize := 6
    #     linecolor := nothing
    #     tipfont := (20, )
    #     markercolor := :white
    #     markerstrokewidth := 0
    #     seriestype := :scatter
    #     # phylo.tree
    #     extrema, [0, 0]
    # end
end