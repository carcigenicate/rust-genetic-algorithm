use rand::{rng, Rng};
use rand::rngs::ThreadRng;
use crate::genetic_algorithm::generate_random_population;

mod genetic_algorithm;

fn main() {
    type GeneT = i32;
    let population_size = 100;

    let fitness_function: genetic_algorithm::FitnessFunction<GeneT> = | genes | {
        let mut even = false;
        let mut fitness: f64 = 0.0;
        for gene in genes.iter() {
            let target_remainder = if even { 0 } else { 1 };
            if gene % 2 == target_remainder {
                fitness += *gene as f64;
            } else {
                fitness -= *gene as f64;
            }
            even = !even;
        }
        fitness
    };

    let gene_producer: genetic_algorithm::GeneProducer<GeneT> = | rand_gen: &mut ThreadRng | {
        rand_gen.random_range(0..=100)
    };

    let tweaks = genetic_algorithm::Tweaks {
        elitism_percent: 0.1,
        population_size,
        mutation_chance: 0.01,
    };

    let mut rand_gen = rng();

    let mut env = genetic_algorithm::Environment {
        tweaks,
        population: generate_random_population(population_size, 10, gene_producer, &mut rand_gen),
        rand_gen,
        fitness_function,
        gene_producer,
        tweak_mutator: None,
    };

    env.advance_while(| i, _ | { i < 1_000 });

    env.sort_population_by_fitness();

    let formatted_results: String = env.population.iter().map(| individual | {
        format!("{}: {:?}", individual.fitness_score.unwrap_or(0.0), individual.genes)
    }).collect::<Vec<String>>().join("\n");

    println!("{}", formatted_results);

}
