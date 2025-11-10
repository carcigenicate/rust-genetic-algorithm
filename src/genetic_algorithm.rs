use rand::{Rng};

use rand::rngs::ThreadRng;

type FitnessScore = f64;
pub type FitnessFunction<G> = fn(&Vec<G>) -> FitnessScore;
pub type GeneProducer<G> = fn(&mut ThreadRng) -> G;
pub type TweakMutator = fn(&mut Tweaks) -> ();

type Population<G> = Vec<Individual<G>>;

#[derive(Clone, Debug)]
pub struct Individual<G: Clone> {
    pub genes: Vec<G>,
    pub fitness_score: Option<FitnessScore>
}

impl<G: Clone> Individual<G> {
    fn new_from(genes: Vec<G>) -> Self {
        Self {
            genes,
            fitness_score: None,
        }
    }

    fn uniform_cross_with(&self, other: &Individual<G>, rand_gen: &mut ThreadRng) -> Self {
        let new_genes: Vec<G> = self.genes.iter().zip(other.genes.iter()).map(| (own_gene, other_gene) | {
            if rand_gen.random_bool(0.5) {
                own_gene.clone()
            } else {
                other_gene.clone()
            }
        }).collect();

        Individual::new_from(new_genes)
    }

    fn new_random(genes_length: usize, gene_producer: GeneProducer<G>, rand_gen: &mut ThreadRng) -> Self {
        let mut genes = Vec::with_capacity(genes_length);
        for _ in 0..genes_length {
            genes.push(gene_producer(rand_gen))
        }

        Individual::new_from(genes)
    }
}

pub struct Tweaks {
    pub population_size: usize,
    // pub n_parents: usize,
    pub elitism_percent: f64,
    pub mutation_chance: f64,
}

pub struct Environment<G: Clone> {
    pub population: Population<G>,
    pub tweaks: Tweaks,
    pub fitness_function: FitnessFunction<G>,
    pub gene_producer: GeneProducer<G>,
    pub tweak_mutator: Option<TweakMutator>,
    pub rand_gen: ThreadRng,
}

impl<G: Clone> Environment<G> {
    fn breed_population(&mut self, children_needed: usize) -> Population<G> {
        const N_PARENTS: usize = 2;

        let total_fitness = self.population.iter().fold(0.0, | running_sum, individual| {
            running_sum + individual.fitness_score.unwrap_or(0.0)  // TODO: Safe default?
        });

        let breeding_pool = self.select_random_mates_roulette(total_fitness, N_PARENTS, children_needed);
        let paired = breeding_pool.chunks_exact(N_PARENTS);

        let mut new_population: Population<G> = Vec::with_capacity(children_needed);
        for pair in paired {
            let mut child = pair[0].uniform_cross_with(&pair[1], &mut self.rand_gen);
            self.mutate_individual(&mut child);
            new_population.push(child);
        }

        new_population
    }

    fn select_random_mates_roulette(&mut self, total_fitness: FitnessScore, n_parents_per_cross: usize, children_needed: usize) -> Vec<Individual<G>> {
        let mut bucketed_population: Vec<(f64, &Individual<G>)> = Vec::with_capacity(self.population.len());
        let mut previous_bucket_max = 0.0;

        for individual in self.population.iter() {
            let bucket_max = previous_bucket_max + (individual.fitness_score.unwrap_or(0.0) / total_fitness);
            previous_bucket_max = bucket_max;
            bucketed_population.push((bucket_max, individual));
        }

        let parents_needed = children_needed * n_parents_per_cross;
        let mut breeding_pool: Vec<Individual<G>> = Vec::with_capacity(parents_needed);
        for _ in 0..parents_needed {
            let random_spin: f64 = self.rand_gen.random();
            for (bucket_max, individual) in bucketed_population.iter() {
                if random_spin < *bucket_max {
                    breeding_pool.push((*individual).clone());
                    break;
                }
            }
        }

        breeding_pool
    }

    fn evaluate_population(&mut self) -> () {
        for i in 0..self.population.len() {
            let individual: &mut Individual<G> = self.population.get_mut(i).unwrap();
            if individual.fitness_score.is_none() {
                individual.fitness_score = Some((self.fitness_function)(&individual.genes))
            }
        }
    }

    pub fn sort_population_by_fitness(&mut self) {
        self.population.sort_by( | individual_one, individual_two | {
            let cmp = individual_one.fitness_score.unwrap_or(0.0).total_cmp(&individual_two.fitness_score.unwrap_or(0.0));
            cmp.reverse()
        });
    }

    fn mutate_individual(&mut self, individual: &mut Individual<G>) -> () {
        let Tweaks { mutation_chance, .. } = self.tweaks;

        for i in 0..individual.genes.len() {
            if self.rand_gen.random_bool(mutation_chance) {
                individual.genes[i] = (self.gene_producer)(&mut self.rand_gen);
            }
        }
    }

    fn advance_generation(&mut self) {
        self.evaluate_population();
        self.sort_population_by_fitness();

        // TODO: Eww
        let n_elite = (self.tweaks.elitism_percent * self.tweaks.population_size as f64).min(self.population.len() as f64).floor() as usize;
        let n_children = self.tweaks.population_size - n_elite;
        let elites = self.population[0..n_elite].to_vec();

        let mut new_children = self.breed_population(n_children);
        new_children.splice(0..0, elites);

        self.population = new_children;
    }

    pub fn advance_while(&mut self, continue_predicate: fn(u64, &Environment<G>) -> bool) {
        for generation in 0u64.. {
            if !(continue_predicate)(generation, &self) {
                break;
            }

            if let Some(tweak_mutator) = self.tweak_mutator.as_mut() {
                tweak_mutator(&mut self.tweaks);
            }

            self.advance_generation();
        }

        self.evaluate_population();
        self.sort_population_by_fitness();
    }
}

pub fn generate_random_population<G: Clone>(population_size: usize, genes_length: usize, gene_producer: GeneProducer<G>, rand_gen: &mut ThreadRng) -> Population<G> {
    let mut population = Vec::with_capacity(population_size);
    for _ in 0..population_size {
        let individual = Individual::new_random(genes_length, gene_producer, rand_gen);
        population.push(individual);
    }

    population
}

