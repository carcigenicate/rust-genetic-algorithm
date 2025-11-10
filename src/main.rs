use rand::{rng, Rng};
use rand::rngs::ThreadRng;

mod genetic_algorithm;
mod n_queens_problem;

fn main() {
    let mut env = n_queens_problem::run_simulation();

    env.sort_population_by_fitness();

    let formatted_results: String = env.population.iter().map(| individual | {
        let mut coords: Vec<(i32, i32)> = individual.genes.chunks_exact(2).map(| pair | (pair[0], pair[1])).collect();
        coords.sort();
        
        format!("{}: {:?}", individual.fitness_score.unwrap_or(0.0), coords)
    }).collect::<Vec<String>>().join("\n");

    println!("{}", formatted_results);

}
