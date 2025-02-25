"""
Parallel processing utilities for genome dereplication.
"""

import logging
import concurrent.futures
from typing import List, Dict, Tuple, Callable, Any, TypeVar, Optional

import tqdm

T = TypeVar("T")
R = TypeVar("R")

logger = logging.getLogger(__name__)


def parallel_map(
    func: Callable[[T], R],
    items: List[T],
    max_workers: int,
    chunksize: int = 1,
    show_progress: bool = False,
    desc: str = None,
) -> List[R]:
    """
    Apply a function to each item in a list in parallel using a process pool.

    Args:
        func: Function to apply to each item
        items: List of items to process
        max_workers: Maximum number of worker processes
        chunksize: Number of items per chunk (for ProcessPoolExecutor)
        show_progress: Whether to show a progress bar
        desc: Description for the progress bar

    Returns:
        List of results
    """
    # Adjust workers if there are fewer items than requested workers
    n_workers = min(max_workers, len(items))

    if n_workers <= 0:
        return []

    results = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Map function to items
        future_to_item = {
            executor.submit(func, item): i for i, item in enumerate(items)
        }

        # Process results as they complete
        if show_progress:
            # Use tqdm for progress reporting
            with tqdm.tqdm(
                total=len(items), desc=desc or "Processing", ncols=80
            ) as pbar:
                for future in concurrent.futures.as_completed(future_to_item):
                    try:
                        result = future.result()
                        results.append((future_to_item[future], result))
                        pbar.update(1)
                    except Exception as e:
                        logger.error(
                            f"Error processing item {future_to_item[future]}: {e}"
                        )
                        results.append((future_to_item[future], None))
                        pbar.update(1)
        else:
            # Process without progress bar
            for future in concurrent.futures.as_completed(future_to_item):
                try:
                    result = future.result()
                    results.append((future_to_item[future], result))
                except Exception as e:
                    logger.error(f"Error processing item {future_to_item[future]}: {e}")
                    results.append((future_to_item[future], None))

    # Sort results by original index
    results.sort()
    return [r for _, r in results]


def parallel_starmap(
    func: Callable[..., R],
    args_list: List[Tuple],
    max_workers: int,
    show_progress: bool = False,
    desc: str = None,
) -> List[R]:
    """
    Apply a function to each tuple of arguments in a list in parallel using a process pool.
    Similar to itertools.starmap but parallel.

    Args:
        func: Function to apply to each tuple of arguments
        args_list: List of argument tuples
        max_workers: Maximum number of worker processes
        show_progress: Whether to show a progress bar
        desc: Description for the progress bar

    Returns:
        List of results
    """

    # Define a wrapper that unpacks arguments
    def wrapper(args_tuple):
        return func(*args_tuple)

    # Use parallel_map with the wrapper
    return parallel_map(
        wrapper, args_list, max_workers, show_progress=show_progress, desc=desc
    )


class ParallelTaskManager:
    """Manager for parallel task execution with progress tracking."""

    def __init__(self, max_workers: int, show_progress: bool = True):
        """
        Initialize parallel task manager.

        Args:
            max_workers: Maximum number of worker processes
            show_progress: Whether to show progress bars
        """
        self.max_workers = max_workers
        self.show_progress = show_progress

    def map(self, func: Callable[[T], R], items: List[T], desc: str = None) -> List[R]:
        """
        Apply a function to each item in a list in parallel.

        Args:
            func: Function to apply to each item
            items: List of items to process
            desc: Description for the progress bar

        Returns:
            List of results
        """
        return parallel_map(
            func, items, self.max_workers, show_progress=self.show_progress, desc=desc
        )

    def starmap(
        self, func: Callable[..., R], args_list: List[Tuple], desc: str = None
    ) -> List[R]:
        """
        Apply a function to each tuple of arguments in a list in parallel.

        Args:
            func: Function to apply to each tuple of arguments
            args_list: List of argument tuples
            desc: Description for the progress bar

        Returns:
            List of results
        """
        return parallel_starmap(
            func,
            args_list,
            self.max_workers,
            show_progress=self.show_progress,
            desc=desc,
        )

    def execute_with_dependencies(
        self,
        tasks: Dict[str, Callable[[], R]],
        dependencies: Dict[str, List[str]],
        desc: str = None,
    ) -> Dict[str, R]:
        """
        Execute tasks in parallel respecting dependencies.

        Args:
            tasks: Dictionary of task ID to task function
            dependencies: Dictionary of task ID to list of dependent task IDs
            desc: Description for the progress bar

        Returns:
            Dictionary of task ID to result
        """
        # Build dependency graph
        dependency_graph = {
            task_id: set(dependencies.get(task_id, [])) for task_id in tasks
        }

        # Track completed tasks and their results
        completed_tasks = set()
        results = {}

        # Initialize progress bar if requested
        total_tasks = len(tasks)
        pbar = None

        if self.show_progress:
            pbar = tqdm.tqdm(
                total=total_tasks, desc=desc or "Executing tasks", ncols=80
            )

        # Process tasks until all are completed
        while len(completed_tasks) < total_tasks:
            # Find tasks that can be executed (all dependencies satisfied)
            ready_tasks = [
                task_id
                for task_id, deps in dependency_graph.items()
                if task_id not in completed_tasks and deps.issubset(completed_tasks)
            ]

            if not ready_tasks:
                if len(completed_tasks) < total_tasks:
                    # This should not happen if dependencies are correctly specified
                    logger.error(
                        "No tasks ready to execute but not all tasks completed. Check for circular dependencies."
                    )
                    break

            # Execute ready tasks in parallel
            task_results = self.map(
                lambda task_id: tasks[task_id](),
                ready_tasks,
                desc=None,  # No nested progress bar
            )

            # Update completed tasks and results
            for task_id, result in zip(ready_tasks, task_results):
                completed_tasks.add(task_id)
                results[task_id] = result

                # Update progress bar
                if pbar:
                    pbar.update(1)

        # Close progress bar
        if pbar:
            pbar.close()

        return results
