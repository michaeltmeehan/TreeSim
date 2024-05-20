function constrain_coalescences!(const_lower::Vector{Float64},
                                 const_upper::Vector{Float64},
                                 const_lineages::Vector{Int64},
                                 const_events::Vector{Int64},
                                 lineages::Vector{Int64},
                                 times::Vector{Float64},
                                 leaves::Vector{Int64},
                                 Nₑ::Float64,
                                 bound::Float64,
                                 norm_tol::Float64)::Float64
    # Total the leaves
    total_leaves = sum(leaves)
    # Likelihood increment
    likelihood = 1.0
    # Initial constraints, bound interval
    events = lineages[2] - lineages[1]
    
    # Counter
    c = 1
    for m in 1:events
        const_lower[c] = bound
        const_upper[c] = times[1]
        const_lineages[c] = lineages[2]
        const_events[c] = events
        c += 1
    end

    for k in 2:length(times)
        events = leaves[k-1] + lineages[k+1] - lineages[k]

        for m in 1:events
            const_lower[c] = times[k-1]
            const_upper[c] = times[k]
            const_lineages[c] = lineages[k+1]
            const_events[c] = events
            c += 1
        end
    end

    # Separate coalescence events
    for c in 1:(total_leaves - 1)
        while (const_events[c] > 1)
            events = const_events[c]
            dt = 0.5 * (const_upper[c] - const_lower[c])
            mid_time = const_upper[c] - dt
            prob_norm = homochronous_probability(const_lineages[c], const_lineages[c] - events, 2.0 * dt, Nₑ)
            sig_loss = significance_loss(const_lineages[c], const_lineages[c] - events, dt, Nₑ)
            if sig_loss > norm_tol
                u = rand()
                events_lhs = 0
                events_rhs = events - events_lhs
                prob_rhs = homochronous_probability(const_lineages[c], const_lineages[c] - events_rhs, dt, Nₑ)
                prob_lhs = homochronous_probability(const_lineages[c] - events_rhs, const_lineages[c] - events, dt, Nₑ)
                sum_prob = (prob_lhs * prob_rhs) / prob_norm
                while u > sum_prob
                    events_lhs += 1
                    events_rhs -= 1
                    prob_rhs = homochronous_probability(const_lineages[c], const_lineages[c] - events_rhs, dt, Nₑ)
                    prob_lhs = homochronous_probability(const_lineages[c] - events_rhs, const_lineages[c] - events, dt, Nₑ)
                    sum_prob += (prob_lhs * prob_rhs) / prob_norm
                end
                likelihood *= (prob_lhs * prob_rhs) / prob_norm
            else
                events_lhs = trunc(Int, events / 2)
                events_rhs = events - events_lhs
                likelihood = 0.0
            end

            for m in 0:(events - 1)
                if m < events_lhs
                    const_upper[c + m] = mid_time
                    const_lineages[c + m] -= events_rhs
                    const_events[c + m] = events_lhs
                else
                    const_lower[c + m] = mid_time
                    const_events[c + m] = events_rhs
                end
            end
        end
    end
    return likelihood
end


constrain_coalescences!(const_lower::Vector{Float64},
                        const_upper::Vector{Float64},
                        const_lineages::Vector{Int64},
                        const_event::Vector{Int64},
                        lineages::Vector{Int64},
                        times::Vector{Float64},
                        leaves::Vector{Int64},
                        Nₑ::Float64,
                        bound::Float64)::Float64 = constrain_coalescences!(const_lower::Vector{Float64},
                                                                           const_upper::Vector{Float64},
                                                                           const_lineages::Vector{Int64},
                                                                           const_event::Vector{Int64},
                                                                           lineages::Vector{Int64},
                                                                           times::Vector{Float64},
                                                                           leaves::Vector{Int64},
                                                                           Nₑ::Float64,
                                                                           bound::Float64,
                                                                           1.0e-10)::Float64




function significance_loss(n_start::Int64, 
                           n_end::Int64, 
                           dt::Float64, 
                           Nₑ::Float64)::Float64
    max_increment = 0.0

    # Check that inputs are valid
    (n_start <= 0 || n_end <= 0 || n_start < n_end || dt < 0 || Nₑ <= 0) && return 1.0
    (n_start == 1 && n_end == 1) && return 1.0

    if n_end == 1
        # Initialise total probability
        prob_total = 0.0
        for k in 2:n_start
            dk = convert(Float64, k)
            # Initialise probability increment
            prob_increment = 1.0
            for l in 2:n_start
                dl = convert(Float64, l)
                if l != k
                    prob_increment *= (dl * (dl - 1.0)) /
                                                (dl * (dl - 1.0) - dk * (dk - 1.0))
                end
            end

            # Calculate exponential
            prob_increment *= 1.0 - exp(- ((dk * (dk - 1.0)) / (2.0 * Nₑ)) * dt)

            if max_increment < abs(prob_increment)
                max_increment = abs(prob_increment)
            end

            # Update total probability
            prob_total += prob_increment
        end
    else
        # Initialise total probability
        prob_total = 0.0
        for k in n_end:n_start
            dk = convert(Float64, k)
            # Initialise probability increment
            prob_increment = (dk * (dk - 1.0)) / (convert(Float64, n_end) * (convert(Float64, n_end) - 1.0))
            for l in n_end:n_start
                dl = convert(Float64, l)
                # Calculate coefficients
                if l != k
                    prob_increment *= (dl * (dl - 1.0)) /
                                                (dl * (dl - 1.0) - dk * (dk - 1.0))
                end
            end

            # Calculate exponential
            prob_increment *= exp(- ((dk * (dk - 1.0)) / (2.0 * Nₑ)) * dt)

            if max_increment < abs(prob_increment)
                max_increment = abs(prob_increment)
            end

            # Update total probability
            prob_total += prob_increment
        end
    end
    return prob_total / max_increment
end