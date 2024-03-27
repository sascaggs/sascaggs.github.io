# plot stuff
library(tidyverse)
library(patchwork)

# network stuff
library(tidygraph)
library(ggraph)
library(igraph)

# link functions 
logit = function(p) {  log( p / (1 - p))  }
inv_logit = function(x) { 1 / (1+ exp(-x) ) }

set.seed(7)

df_theme =   
    theme(plot.background = element_rect(fill='white'), 
          panel.background = element_rect(fill='black'), 
          panel.grid = element_blank())

plot_theme =   
    theme(plot.background = element_rect(fill='white'), 
          panel.background = element_rect(fill='white', color='black'), 
          panel.grid = element_blank())

# testing ----
N = 25

# create species attributes 
d = data.frame(
    species = 1:N, 
    troph = sample( c('carn','omni','herb','auto'), size = N, replace = TRUE, prob = c(0.3,0.05,0.25,0.4) ), 
    body_size = rgamma(N, shape = 2))
d

# create and tidy trophic links 
g = expand.grid(1:N, 1:N)
g = merge(g, d, by.x = 'Var1', by.y = 'species')
g = merge(g, d, by.x = 'Var2', by.y = 'species')

# merging leads to function orders of columns, rename carefully 
colnames(g) = c('prey','pred','troph_pred','body_size_pred', 'troph_prey', 'body_size_prey')

# tidy it up 
g = g %>% 
    arrange(pred,prey) %>%
    mutate(dyad = 1:nrow(g)) %>%
    select(pred, prey, dyad, troph_pred, troph_prey, body_size_pred, body_size_prey) 
g

# calculate body size differences 
g$body_size_diff = (g$body_size_pred - g$body_size_prey)+0.1 

# carnivores 
g$carnivory = ifelse(g$troph_pred == 'carn' & !g$troph_prey == 'auto', 1, 0)

# omnivores 
g$omnivory = ifelse(g$troph_pred == 'omni' , 1, 0)

# herbivory 
g$herbivory = ifelse(g$troph_pred == 'herb' & g$troph_prey == 'auto'  , 1, 0)


bC = 1.5
bO = 1.5
bH = 1.5
bS = -1.5

inv_logit( bC*g$herbivory + bO*g$omnivory + bH*g$herbivory )




# functionalized food web simulation ----

# key
# - N: number of species
# - seed: reproducibility
# - weights: samples weights for carnivore, ominovore, herbivore, and autotrophs

sim_FW = function( N , 
                   seed = 777 ,  
                   weights = c(0.3, 0.05, 0.25, 0.4) , 
                   shape = 2 , 
                   model = 'allometric' , 
                   base_rate = -2 , 
                   alpha , 
                   bC , 
                   bO ,  
                   bH ) {
    
    set.seed(seed)
    # trophic levels 
    trophs = c('carn','omni','herb','auto')
    
    if ( length(weights) != length(trophs)) { print('Error: length of probability vector must equal 4.') } else {
        
        # create species attributes 
        d = data.frame(
            species = 1:N, 
            troph = sample( trophs, size = N, replace = TRUE, prob = weights ), 
            body_size = rgamma(N, shape = shape)
        ) 
    }
    
    # create and tidy trophic links 
    g = expand.grid(1:N, 1:N)
    g = merge(g, d, by.x = 'Var1', by.y = 'species')
    g = merge(g, d, by.x = 'Var2', by.y = 'species')
    
    # merging leads to function orders of columns, rename carefully 
    colnames(g) = c('prey','pred','troph_pred','body_size_pred', 'troph_prey', 'body_size_prey')
    
    # tidy it up 
    g = g %>% 
        arrange(pred,prey) %>%
        mutate(dyad = 1:nrow(g)) %>%
        select(pred, prey, dyad, troph_pred, troph_prey, body_size_pred, body_size_prey)
    
    # calculate body size differences 
    g$body_size_diff = (g$body_size_pred - g$body_size_prey )  + 0.1
    
    # carnivores 
    g$carnivory = ifelse(g$troph_pred == 'carn' & !g$troph_prey == 'auto', 1, 0)
    
    # omnivores 
    g$omnivory = ifelse(g$troph_pred == 'omni' , 1, 0)
    
    # herbivory 
    g$herbivory = ifelse(g$troph_pred == 'herb' & g$troph_prey == 'auto'  , 1, 0)
    
    # compute probabilities 
    # impossible trophic interactions are set to 0 (e.g, a carnivore browsing on grass)
    
    if( model == 'allometric') { 
        
        D = g$body_size_diff
        g$p = inv_logit( base_rate + sign(D)*abs(D)^alpha  )
        g$consumption = rbinom( nrow(g), size=1, prob = g$p ) 
        g[ g$troph_pred == 'auto' , 'consumption' ] = 0
        g[ g$troph_pred == 'carn' & g$troph_prey == 'auto' , "consumption"] = 0
        
        
    } else if( model == 'dietary' ) {
        
        g$p = inv_logit( base_rate + bC*g$carnivory + bO*g$omnivory + bH*g$herbivory )
        g$consumption = rbinom( nrow(g), size=1, prob = g$p ) 
        g[ g$troph_pred == 'auto' , 'consumption' ] = 0
        g[ g$troph_pred == 'carn' & g$troph_prey == 'auto' , "consumption"] = 0
        
        
    } else if( model == 'combined' ) {
        
        D = g$body_size_diff
        g$p = inv_logit( base_rate + sign(D)*abs(D)^alpha + bC*g$carnivory + bO*g$omnivory + bH*g$herbivory )
        g$consumption = rbinom( nrow(g), size=1, prob = g$p ) 
        g[ g$troph_pred == 'auto' , 'consumption' ] = 0
        g[ g$troph_pred == 'carn' & g$troph_prey == 'auto' , "consumption"] = 0
        
        
    }
    
    el = g[ g$consumption == 1 ,c('pred','prey','consumption','p','body_size_diff') ] 
    ig = graph.data.frame( el, directed = T, vertices = d)

    d$generality = graph.strength(ig, mode = 'out')
    d$vulnerability = graph.strength(ig, mode = 'in')
    
    d$outcore = graph.coreness(ig, mode = 'out')
    d$incore = graph.coreness(ig, mode = 'in')
    
    list(d,g,ig)

}





# test plots ----
sim1 = sim_FW(N = 100, weights = c(0.1,0.1,0.1,0.7), seed = NULL,
              alpha = 1.13, base_rate = -3, model = 'allometric' )


sim1[[2]] %>% 
    filter(troph_pred == 'carn' & troph_prey == 'auto') %>%
    select(consumption) %>%
    summarise(s = unique(consumption))


el = sim1[[2]] %>% 
    select(pred, prey, body_size_diff, consumption ) %>%
    filter(consumption == 1) 

ig = graph.data.frame(el, directed = F, vertices = sim1[[1]])


ig %>%
    as_tbl_graph() %>%
    ggraph('mds') + 
    geom_edge_link(color='#00000066', width=0.681) + 
    geom_node_point(aes(fill=troph, size=body_size), shape=21, stroke = 1.618) + 
    theme(plot.background = element_rect(fill='white'), 
          panel.background = element_rect(fill='white')) + 
    scale_fill_manual(values = c('#3300ff','#ff0066','#99ff99','#ffaa00'))





cl = cluster_leiden(ig, objective_function = 'modularity')
cl$membership

modularity(ig, membership = cl$membership, directed = F)


#description 

tibble(generality = graph.strength(ig, mode = 'out'), 
       troph = V(ig)$troph, 
       body_size = V(ig)$body_size ) %>%
    filter(!troph == 'a') %>%
    ggplot() + geom_point(aes(x=body_size, y=generality))


######### ALLOMETRIC MODEL ----
N_iter = 100
sims = replicate( n=N_iter, simplify = 'array',
                  expr =  sim_FW(N=50, weights = c(0.1,0.1,0.1,0.7), seed = NULL,
                                 alpha = 1.13, base_rate = -4, model = 'allometric' ) )

#sims[1,][[1]][,"iter"] = 1

# assign iteration ids 
for(i in 1:N_iter) { for(j in 1:2)  { sims[j,][[i]][, "iter" ] = i } }

# combine tables
sim_df = bind_rows( t(sims)[,1] )
sim_dy = bind_rows( t(sims)[,2] )

# summarise graphs 
ig_sims = sims[3,]

sims_sum = data.frame(
    iter = 1:N_iter,
    model = 'allometric',
    vertices = unlist(lapply(ig_sims, vcount)), 
    edges = unlist(lapply(ig_sims, ecount)), 
    density = unlist(lapply(ig_sims, graph.density)), 
    apl = unlist(lapply(ig_sims, average.path.length)), 
    reciprocity = unlist(lapply(ig_sims, reciprocity)), 
    diameter = unlist(lapply(ig_sims, igraph:::diameter)), 
    transitivity = unlist(lapply(ig_sims, igraph::transitivity, type = 'average'))
    
)

sims_sum %>%
    select(density:transitivity) %>%
    gather() %>%
    ggplot(aes(x=value)) + geom_density(color='white', fill='white', alpha=0.5) +
    facet_wrap(~key, scale='free') + 
    df_theme



sim_df %>%
    #filter(!troph == 'auto') %>%
    ggplot(aes(x=body_size, y=generality)) + 
    df_theme + 
    geom_point(aes(color=troph), size=2, stroke=2, alpha=0.25, shape=21 )  +
    scale_color_manual(values = c('#3300ff','#ff0066','#99ff99','#ffaa00')) + 
    ylim(c(0,45)) + 

sim_df %>%
    #filter(!troph == 'auto') %>%
    ggplot(aes(x=body_size, y=vulnerability)) +
    df_theme + 
    geom_point(aes(color=troph), size=2, stroke=2, alpha=0.25, shape=21 )  +
    scale_color_manual(values = c('#3300ff','#ff0066','#99ff99','#ffaa00')) + 
    ylim(c(0,45)) 



######### DIETARY MODEL ----
N_iter = 100
diet_sims = replicate( n=N_iter, simplify = 'array',
                  expr =  sim_FW(N=50, model = 'dietary', weights = c(0.1,0.1,0.1,0.7), seed = NULL,
                                 base_rate = -4, bC=0.25, bO=0.5, bH=0.5 ) )

#diet_sims[1,][[1]][,"iter"] = 1

# assign iteration ids 
for(i in 1:N_iter) { for(j in 1:2)  { diet_sims[j,][[i]][, "iter" ] = i } }

diet_df = bind_rows( t(diet_sims)[,1] )
diet_dy = bind_rows( t(diet_sims)[,2] )

ig_diet = diet_sims[3,]


diet_sum = data.frame(
    iter = 1:N_iter,
    model = 'dietary',
    vertices = unlist(lapply(ig_diet, vcount)), 
    edges = unlist(lapply(ig_diet, ecount)), 
    density = unlist(lapply(ig_diet, graph.density)), 
    apl = unlist(lapply(ig_diet, average.path.length)), 
    reciprocity = unlist(lapply(ig_diet, reciprocity)), 
    diameter = unlist(lapply(ig_diet, igraph:::diameter)), 
    transitivity = unlist(lapply(ig_diet, igraph::transitivity, type = 'average'))
    
)

diet_sum %>%
    select(density:transitivity) %>%
    gather() %>%
    ggplot(aes(x=value)) + geom_density(color='white', fill='white', alpha=0.5) +
    facet_wrap(~key, scale='free') + 
    df_theme


diet_df %>%
    filter(!troph == 'auto') %>%
    ggplot(aes(x=body_size, y=generality)) + 
    geom_point(aes(color=troph), size=2, stroke=2, alpha=0.2, shape=21 ) + 
    df_theme + 
    geom_point(aes(color=troph), size=2, stroke=2, alpha=0.25, shape=21 )  +
    scale_color_manual(values = c('#3300ff','#ff0066','#99ff99','#ffaa00')) + 
    ylim(c(0,12)) +
    
    diet_df %>%
    #filter(!troph == 'auto') %>%
    ggplot(aes(x=body_size, y=vulnerability)) + 
    df_theme + 
    geom_point(aes(color=troph), size=2, stroke=2, alpha=0.25, shape=21 )  +
    scale_color_manual(values = c('#3300ff','#ff0066','#99ff99','#ffaa00')) + 
    ylim(c(0,12)) 

######### COMBINED MODEL ----
N_iter = 100
comb_sims = replicate( n=N_iter, simplify = 'array',
                       expr =  sim_FW(N=50, model = 'combined', weights = c(0.1,0.1,0.1,0.7), seed = NULL,
                                      base_rate = -4, alpha=2,  bC=0.25, bO=0.5, bH=0.5 ) )

#diet_sims[1,][[1]][,"iter"] = 1

# assign iteration ids 
for(i in 1:N_iter) { for(j in 1:2)  { diet_sims[j,][[i]][, "iter" ] = i } }

comb_df = bind_rows( t(diet_sims)[,1] )
comb_dy = bind_rows( t(diet_sims)[,2] )

ig_comb = comb_sims[3,]


diet_sum = data.frame(
    iter = 1:N_iter,
    model = 'dietary',
    vertices = unlist(lapply(ig_comb, vcount)), 
    edges = unlist(lapply(ig_comb, ecount)), 
    density = unlist(lapply(ig_comb, graph.density)), 
    apl = unlist(lapply(ig_comb, average.path.length)), 
    reciprocity = unlist(lapply(ig_comb, reciprocity)), 
    diameter = unlist(lapply(ig_comb, igraph:::diameter)), 
    transitivity = unlist(lapply(ig_comb, igraph::transitivity, type = 'average'))
    
)

diet_sum %>%
    select(density:transitivity) %>%
    gather() %>%
    ggplot(aes(x=value)) + geom_density(color='white', fill='white', alpha=0.5) +
    facet_wrap(~key, scale='free') + 
    df_theme


comb_df %>%
    filter(!troph == 'auto') %>%
    ggplot(aes(x=body_size, y=generality)) + 
    df_theme + 
    geom_point(aes(color=troph), size=2, stroke=2, alpha=0.25, shape=21 )  +
    scale_color_manual(values = c('#3300ff','#ff0066','#99ff99','#ffaa00')) + 
    ylim(c(0,12)) + 
    
    comb_df %>%
    #filter(!troph == 'auto') %>%
    ggplot(aes(x=body_size, y=vulnerability)) + 
    df_theme + 
    geom_point(aes(color=troph), size=2, stroke=2, alpha=0.25, shape=21 )  +
    scale_color_manual(values = c('#3300ff','#ff0066','#99ff99','#ffaa00')) + 
    ylim(c(0,12)) 



#plot( seq(-2,2, length=9), log(exp(seq(-2,2, length=9))^1), xlim = c(-2,2), ylim = c(-3,3), type='l' )
#lines( seq(-2,2, length=9), log(exp(seq(-2,2, length=9))^1.5), xlim = c(-2,2), ylim = c(-3,3), type='l', col='red')
#lines( seq(-2,2, length=9), log(exp(seq(-2,2, length=9))^0.5), xlim = c(-2,2), ylim = c(-3,3), type='l', col='blue')
#
#
#plot( seq(-2,2, length=9), exp(seq(-2,2, length=9))^1, xlim = c(-2,2), ylim = c(0,3), type='l' )
#lines( seq(-2,2, length=9), exp(seq(-2,2, length=9))^1.5, xlim = c(-2,2), ylim = c(0,3), type='l', col='red')
#lines( seq(-2,2, length=9), exp(seq(-2,2, length=9))^0.5, xlim = c(-2,2), ylim = c(0,3), type='l', col='blue')

