#include <stdio.h>
#include <complex.h>
#include <math.h>

/* * A PEDANTIC NOTE ON HILBERT SPACES:
 * We define the qubit not merely as a data type, but as a normalized vector 
 * residing within a complex projective space CP^1.
 * * The 'flex' here is recognizing that standard C pointers are insufficient 
 * to represent 'spooky action at a distance' without explicit shared state 
 * simulation (a tensor product representation).
 */

#define ONE_OVER_ROOT_2 (1.0 / sqrt(2.0))

// Defining the Basis States
typedef struct {
    double complex alpha; // Amplitude for |0>
    double complex beta;  // Amplitude for |1>
} Qubit;

// A Two-Qubit System (The Tensor Product Space H_A (x) H_B)
typedef struct {
    double complex c00; // |00>
    double complex c01; // |01>
    double complex c10; // |10>
    double complex c11; // |11>
} EntangledSystem;

// The Hadamard Gate (H): Maps |0> -> (|0> + |1>)/sqrt(2)
// This creates the initial superposition required for the flex.
Qubit apply_hadamard(Qubit q) {
    return (Qubit){
        .alpha = ONE_OVER_ROOT_2 * (q.alpha + q.beta),
        .beta  = ONE_OVER_ROOT_2 * (q.alpha - q.beta)
    };
}

// The CNOT Gate: Flips target if control is |1>.
// This is the "Parenthetical" operation where independent systems lose individual identity.
EntangledSystem entangle(Qubit control, Qubit target) {
    // Mathematically: CNOT |x, y> = |x, x XOR y>
    // We compute the tensor product (control (x) target) first, then permit the interaction.
    
    EntangledSystem system;
    
    // Applying CNOT logic to the amplitudes of the tensor product:
    // |00> remains |00> (control is 0)
    system.c00 = control.alpha * target.alpha;
    
    // |01> remains |01> (control is 0)
    system.c01 = control.alpha * target.beta;
    
    // |10> becomes |11> (control is 1, target flips 0->1)
    // Note: We swap the coefficients for the 'control=1' subspace.
    system.c11 = control.beta * target.alpha; 
    
    // |11> becomes |10> (control is 1, target flips 1->0)
    system.c10 = control.beta * target.beta;

    return system;
}

int main() {
    /* * STEP 1: INITIALIZATION
     * We begin with a separable state |Psi> = |0> (x) |0>.
     * This is the mundane, classical reality we are attempting to escape.
     */
    Qubit qA = { .alpha = 1.0 + 0.0*I, .beta = 0.0 + 0.0*I }; // |0>
    Qubit qB = { .alpha = 1.0 + 0.0*I, .beta = 0.0 + 0.0*I }; // |0>

    /*
     * STEP 2: SUPERPOSITION
     * Apply Hadamard to qA. 
     * qA becomes |+> = 1/sqrt(2)(|0> + |1>).
     * The system is now (1/sqrt(2)(|0> + |1>)) (x) |0>.
     */
    qA = apply_hadamard(qA);

    /* * STEP 3: ENTANGLEMENT (The Flex)
     * We apply the CNOT gate using qA as control and qB as target.
     * The state evolves from a separable product state to the Bell State |Phi+>:
     * |Phi+> = 1/sqrt(2) * (|00> + |11>)
     * * At this point, writing 'qA' or 'qB' individually is mathematically incoherent.
     * They no longer possess individual state vectors; they are parenthetically bound.
     */
    EntangledSystem bell_state = entangle(qA, qB);

    // Outputting the Probability Amplitudes of the Bell State
    printf("--- The Bell State |Phi+> ---\n");
    printf("|00> Amplitude: %.3f + %.3fi (Prob: %.1f%%)\n", 
           creal(bell_state.c00), cimag(bell_state.c00), cabs(bell_state.c00)*cabs(bell_state.c00)*100);
           
    printf("|01> Amplitude: %.3f + %.3fi (Prob: %.1f%%)\n", 
           creal(bell_state.c01), cimag(bell_state.c01), cabs(bell_state.c01)*cabs(bell_state.c01)*100);

    printf("|10> Amplitude: %.3f + %.3fi (Prob: %.1f%%)\n", 
           creal(bell_state.c10), cimag(bell_state.c10), cabs(bell_state.c10)*cabs(bell_state.c10)*100);

    printf("|11> Amplitude: %.3f + %.3fi (Prob: %.1f%%)\n", 
           creal(bell_state.c11), cimag(bell_state.c11), cabs(bell_state.c11)*cabs(bell_state.c11)*100);
           
    printf("\nStatus: Local Realism Violated.\n");

    return 0;
}
