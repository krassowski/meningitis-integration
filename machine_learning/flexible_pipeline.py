from pandas import Series
from sklearn.pipeline import Pipeline, make_pipeline


def is_instance_by_name(child, parent):
    """Like isinstance, but works with autoreloader.

    Unfortunately, is sensitive to location (will yield false positives).
    """
    return parent.__name__ in {x.__name__ for x in child.__class__.mro()}


class FlexiblePipeline(Pipeline):
    """Pipeline that can be easily modified in the runtime,

    so that the heavy lifting steps (data munging etc),
    can be executed first and excluded from further steps,
    which is useful for cross-validation.

    A simple alternative is to not include data munging like
    steps in the pipeline in the first place, but then it may
    get difficult to keep track of everything - this wrapper
    allows to keep the full pipeline in a single place.
    """

    def apply_steps(self, kind, data):
        for step, transformer in self.steps:
            if is_instance_by_name(transformer, kind):
                data = transformer.fit_transform(data)
        return data

    def exclude(self, to_exclude):
        return make_flexible_pipeline(
            *[
                transformer
                for step, transformer in self.steps
                # again, a workaround for autoreloader
                # if type(transformer) not in to_exclude
                if transformer.__class__.__name__ not in {
                    t.__name__
                    for t in to_exclude
                }
            ]
        )

    def get_coefficients(self, index: Series):
        estimator = self._final_estimator
        return Series(estimator.coef_[0], index=index)


def make_flexible_pipeline(*args, **kwargs):
    pipeline = make_pipeline(*args, **kwargs)
    pipeline.__class__ = FlexiblePipeline
    return pipeline
