use std::collections::HashSet;

use rand::{rng, Rng};
use rand::rngs::ThreadRng;
use crate::genetic_algorithm;
use crate::genetic_algorithm::generate_random_population;

type GeneT = i32;
type Coord = (GeneT, GeneT);
type IndexedCoord = (usize, Coord);

fn count_overlapping(side_length: i32, current_indexed_coord: IndexedCoord, all_coords: &Vec<IndexedCoord>) -> usize {
    let (current_i, (current_x, current_y)) = current_indexed_coord;
    // A way to generate this once? How to deal with the current value?
    let relevant_coords = all_coords.iter().filter_map(| (i, coord ) | if *i == current_i { None } else { Some(coord) } );
    let other_coord_set = HashSet::<_>::from_iter(relevant_coords);
    let mut count = 0;
    for n in 0i32..side_length {
        if other_coord_set.contains(&(n, current_y)) {
            count += 1;
        }

        if other_coord_set.contains(&(current_x, n)) {
            count += 1;
        }

        if other_coord_set.contains(&(current_x + n, current_y + n))
           || other_coord_set.contains(&(current_x + n, current_y - n))
           || other_coord_set.contains(&(current_x - n, current_y + n))
           || other_coord_set.contains(&(current_x - n, current_y - n)
        ) {
            count += 1;
        }
    }

    count
}

fn fitness_function(genes: &Vec<GeneT>) -> f64 {
    let coords: Vec<IndexedCoord> = genes.chunks_exact(2)
        .enumerate()
        .map(|(i, pair)| (i, (pair[0], pair[1])))
        .collect();

    let mut total_overlapping = 0;
    for current_indexed_coord in coords.iter() {
        total_overlapping += count_overlapping(8, *current_indexed_coord, &coords);
    }

    -(total_overlapping as f64)
}

pub fn run_simulation() -> crate::genetic_algorithm::Environment<GeneT> {
    let population_size = 10;

    let gene_producer: genetic_algorithm::GeneProducer<GeneT> = | rand_gen: &mut ThreadRng | {
        rand_gen.random_range(0..8)
    };

    let tweaks = genetic_algorithm::Tweaks {
        elitism_percent: 5.0 / population_size as f64,
        population_size,
        mutation_chance: 0.1,
    };

    let mut rand_gen = rng();

    let mut env = genetic_algorithm::Environment {
        tweaks,
        population: generate_random_population(population_size, 16, gene_producer, &mut rand_gen),
        rand_gen,
        fitness_function,
        gene_producer,
        tweak_mutator: None,
    };

    env.advance_while(| generation, env | {
        let best_fitness = env.population[0].fitness_score.unwrap_or(f64::NEG_INFINITY);

        if generation % 10_000 == 0 {
            println!("Generation: {}, Best Fitness: {}, Best Genes: {:?}", generation, best_fitness, env.population[0].genes);
        }

        let found_solution = best_fitness >= 0.0;
        let attempts_exceeded = generation > 2_000_000;

        !attempts_exceeded && !found_solution
    });

    env
}

// Found:
//  - [(1, 1), (5, 7), (7, 0), (2, 3), (3, 6), (6, 5), (0, 4), (4, 2)]
//  - [(2, 6), (0, 5), (4, 3), (7, 4), (3, 1), (5, 7), (1, 2), (6, 0)]
//  - [(2, 7), (1, 0), (3, 5), (6, 1), (5, 6), (4, 2), (0, 4), (7, 3)]
//  - [(4, 6), (3, 0), (2, 3), (1, 7), (6, 5), (5, 1), (7, 2), (0, 4)]
//  - [(0, 3), (2, 6), (4, 5), (7, 0), (3, 2), (6, 4), (1, 1), (5, 7)]
//  - [(0, 3), (1, 5), (2, 7), (3, 1), (4, 6), (5, 0), (6, 2), (7, 4)]  // Symmetrical solution used as the Wikipedia example
//  - [(0, 1), (1, 4), (2, 6), (3, 3), (4, 0), (5, 7), (6, 5), (7, 2)]
//  - [(0, 2), (1, 5), (2, 3), (3, 0), (4, 7), (5, 4), (6, 6), (7, 1)]
//  - [(0, 6), (1, 4), (2, 2), (3, 0), (4, 5), (5, 7), (6, 1), (7, 3)]
//  - [(0, 5), (1, 3), (2, 6), (3, 0), (4, 7), (5, 1), (6, 4), (7, 2)]